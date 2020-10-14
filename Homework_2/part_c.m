%% AuE 893: Robust Predictive Conrtol HW02, Question 2.18
% Author: Vipul Kumbhar
% ref:

clc;clear;close all;

% All important variables
time_step      = 0.5;                % [sec]

% Set points for MPC
y1_set         = 1;                  % Unit setpoint change in first output 
y2_set         = 0;                  % No setpoint change in first output
mpc_time       = 40;                 % [sec]
Tph            = 10;                 % [sec] prediction horizon

%% Given system model

S_tf    = tf('s');
s_a     = 2/((10*S_tf)+1);
s_b     = 1/(S_tf+1);
s_c     = 2/(S_tf+1);
s_d     = -4/(S_tf+1);
H       = [s_a s_b;s_c s_d];

tf_continous   = ss(H);

% Conti. to descrete function
tf_descrete    = c2d(tf_continous,time_step);

% Print descrete function model
% fprintf("The given model is \n")
sys            = ss(tf_descrete)

% Compare continuous and descrete function response
% step(tf_continous,'r:',tf_descrete,'--');

% State space model
A              = sys.A;
B              = sys.B;
C              = sys.C;
D              = sys.D;

%%
% optimized set points X1,X2,X3
Y_setpoint     = [y1_set;y2_set];
I              = eye(3);
Aeq            = [I-A -B; C D];           % Augmented matrix
beq            = [zeros(3,1); Y_setpoint];% b, in Ax=b
X_guess_sp     = 0.5*ones(5,1);           % initial values,x0, u0

% optimization (useful when output is not feasible/gives lowest cost value) 
objective      = @objc_218;               % objective function
Xopt           = fmincon(@(x)objective(x,C,Y_setpoint),X_guess_sp,[],[],Aeq,beq,[],[],[]);

% Setpoints triple
ys_feasible_3  = C*Xopt(1:3,:);
X_setpoint     = Xopt(1:3,:);
Usp            = Xopt(4:5,:);

if ys_feasible_3 == Y_setpoint
    fprintf('Setpoint is feasible')
end

%%
% MPC prediction 
Hp             = Tph / time_step;         % Piecewise intervals

% Initialization 
t_initial      = 20;
t_end          = t_initial +Tph;          % [sec]

% plant model simulation for first t_initial seconds
u_sim2         = Usp.*ones(2,((t_end/time_step +1)));
u_sim          = zeros(2,(t_end/time_step +1));
t_sim          = 0:time_step:t_end;
x_ini          = [0 0 0];

% ODE45 equivalent, system response
[y,t,x]        = lsim(sys,u_sim,t_sim,x_ini);

% Plot system response for initial simulation
figure('Name','System Response')
subplot(4,1,1)
plot(t,y)
title('Output vs time plot')
legend('Y1','Y2')

subplot(4,1,2)
plot(t,x)
title('States vs time plot')
legend('X1','X2','X3')

subplot(4,1,3)
plot(t,u_sim)
title('Input vs time plot')
legend('U1','U2')
xlabel('Time (sec)')

[y2,t2,x2] = lsim(sys,u_sim2,t_sim,x_ini);
subplot(4,1,4)
plot(t2,y2)
title('Outputs for given setpoint plot')
legend('Y1','Y2')
xlabel('Time (sec)')
axis([0 t_end -1 y1_set+2])

close all;

%%
% Separing the outputs for states x from simulation outputs
% x - system state obtained from lsim
% y - system output obtained from lsim

% after initial run of system for t_initial time
Y_ini           = (y(end - Hp,:))';
X_ini           = (x(end - Hp,:))';

% initial guesses for the optimization variables for the first prediction
% horizon t_initial to (t_initial+prediction horizon)
x_0             = x((end-Hp+1):end,:)'; % Initial guess for state
u_0             = Usp.*ones(2,Hp+1);    % Initial guess for input

%u_0             = u_sim(2,end - Hp:end);% Initial guess for input
slack_0         = 2;                    % slack variable is added to the ...
                        ...optimization variable for feasibility of the NLP
                            
%initial guess states of all variables [x u slack]                           
x0             = [x_0(1,:) x_0(2,:) x_0(3,:) u_0(1,:) u_0(1,:) slack_0];
time           = 0;                    % MPC start time

%% MPC - NLC loop from t_initial to (t_initial+mpc_time) with time_step
mpc_time_steps     = (mpc_time/time_step);

% reset usim and tsim to begining of mpc control time
u_sim = u_sim(:,1:end-Hp);
t_sim = t_sim(:,1:end-Hp);

U_ini = u_sim(:,end-Hp);

for i=1:mpc_time_steps
    tic
    
    % Linear inequality constrains for x
    ...ref: https://www.mathworks.com/help/optim/ug/linear-constraints.html
        
    % None for a part
    % |u(k)|<1 for b part
    % |delta(U)/delta(t)| < 0.1 for c part
    
    % for part a
    Aeq     = [];
    beq     = [];
    
    % Linear equlaity constrains for x
    Am      = [];
    bm      = [];
    
    % lower and upper bounds for states, u and slack(v)
    lb_x1           = -inf*ones(1,Hp);
    lb_x2           = -inf*ones(1,Hp);
    lb_x3           = -inf*ones(1,Hp);
    lb_u1           = -1*ones(1,Hp+1);
    lb_u2           = -1*ones(1,Hp+1);
    lb_v            = 0;
    
    ub_x1           = inf*ones(1,Hp);
    ub_x2           = inf*ones(1,Hp);
    ub_x3           = inf*ones(1,Hp);
    ub_u1           = 1*ones(1,Hp+1);
    ub_u2           = 1*ones(1,Hp+1);
    ub_v            = inf;
    
    % Upper and lower bounds for all the variable -combined
    lb              = [lb_x1 lb_x2 lb_x3 lb_u1 lb_u2 lb_v];
    ub              = [ub_x1 ub_x2 ub_x3 ub_u1 ub_u2 ub_v];
 
    % fmincon optimization options
    options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
    nonlcon         = @nonlin_c;
    call_fun        = @obj_fun;

    % fmincon NLP solver
    % opt_var - optimized variable output, fval - function value
    
    [opt_var, fval] = fmincon(@(opt)call_fun(opt,time_step,X_setpoint(1),X_setpoint(2),X_setpoint(3),Usp,C,Hp),...
        x0,Am,bm,Aeq,beq,lb,ub,@(opt)nonlcon(opt,X_ini,U_ini,A,B,time_step,Hp),options);
    
    % Initial guess for optimization varibles for the next prediction
    % horizon
    x0              = opt_var;
    
    % first time step control input from fmincon optimization
    no_of_states    = 3;
    U_1             = opt_var((no_of_states*Hp+1));
    U_2             = opt_var(((no_of_states+1)*Hp+2));
    U               = [U_1;U_2];
    
    % Simulate the plant model from new U's to get system states at time
    % t_initial + i*time_step
    % also adds new input to u_sim , used for plots 
    u_sim           = cat(2,u_sim,U);
    t_sim           = cat(2,t_sim,(time_step+t_sim(end)));
    [y,t,x]         = lsim(sys,u_sim,t_sim,[0;0;0]);
    
    % Storing output states for next prediction horizon             
    X_ini           = x(end,:)';U_ini = U;
    
    % Shift time step
    time            = time + time_step;
        
    toc
end
  
%% Plot 
figure('Name','MPC Optimization plots')
subplot(3,1,1)
stairs(t,u_sim(1,:),'LineWidth',2);hold on
stairs(t,u_sim(2,:),'LineWidth',2)
title('Input vs Time plot');legend('U1','U2')
axis([0 t(end) min(min(u_sim))-0.2 max(max(u_sim))+0.2])

subplot(3,1,2)
y_set_plot1 = y1_set*ones(length(t))';
plot(t,y,t,y_set_plot1,':','LineWidth',2)
title('Output states vs time plot')
axis([0 t(end) -1 y1_set+2]);legend('Y1','Y2','Y1 setpoint')

subplot(3,1,3)
x_set_plot1 = X_setpoint(1)*ones(length(t))';
x_set_plot2 = X_setpoint(2)*ones(length(t))';
x_set_plot3 = X_setpoint(3)*ones(length(t))';
plot(t,x,t,x_set_plot1,':',t,x_set_plot2,':',t,x_set_plot3,':','LineWidth',2)
title('States vs time plot')
legend('X1','X2','X3','X1 setpoint','X2 setpoint','X3 setpoint')
axis([0 t(end) -1 max(X_setpoint)+2])
xlabel('Time [Seconds]')
  
%% The End