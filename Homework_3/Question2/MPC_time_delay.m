%% AuE 893: Robust Predictive Conrtol HW03, Excercise 1.61
% Author: Vipul Kumbhar
% ref   : AuE 893 HW2 solutions

clear;clc;close all;tic

% Important parameters
del_t     = 1;                   % Sample time for discretization
t_0       = 0;                   % MPC start time
t_sim     = 80;                  % Total simulation time
N         = 15;                  % Number of intervals
Tp        = del_t*N;             % Prediction horizon length
mpc_steps = (t_sim - t_0)/del_t; % number of steps

y1_sp1    = 1;                   % Y setpoint from t = 0 to t = step_change
step_change= 20;                 % y1 set point step change after 20 steps
y1_sp2    = 1;                   % Y setpoint after t = step_change

% To find the state space of the transfer function
s         = tf('s');
theta     = 5;
k         = 1;
tou       = 1;

% Ref:https://www.mathworks.com/help/mpc/ug/plant-models-with-delays.html
G_s       = tf([k],[tou 1],'IOdelay',theta,'TimeUnit','seconds')
sysc      = ss(G_s);

% Continuous to Discrete State-space form
sysd      = c2d(sysc,del_t,'zoh')
A         = sysd.A;
B         = sysd.B;
C         = sysd.C;
D         = sysd.D;

% Size of the state matrices
n         = size(A,1);
m         = size(B,2);
p         = size(C,1);
na        = n+m;

% Augmented state, input and output matrices
Aa        = blkdiag(A,zeros(m));
Ba        = [B; eye(m)];
Ca        = [C zeros(p,m)];

% Disturbance matrix
nd            = p;          % Number of disturbance = number of outputs
C_disturbance = eye(nd);
B_disturbance = zeros(n,nd);

Cd        = C_disturbance;
% augmented disturbance matrix (B_disturbance = )
Bd        = [B_disturbance;zeros(1,1)];

%
% Detectability check
if (rank([eye(na)-Aa -Bd; Ca Cd]) == na+nd)
    disp('The augmented system is detectable');
else
    disp('The augmented system is not detectable');
end

%%
% State matrices augmented with disturbance
A_tilda   = [Aa Bd; zeros(nd,n+m) eye(nd)];
B_tilda   = [Ba; zeros(nd,m)];
C_tilda   = [Ca Cd];

G_d       = blkdiag(zeros(na),eye(nd));
H_d       = zeros(p,na+nd);

%% System estimate
kalman_sys= ss(A_tilda,[B_tilda G_d],C_tilda,[D H_d],del_t);

% Process and measurement noise covariances
Qw        = blkdiag(zeros(na),0.005*eye(nd));
Rv        = 0.01*eye(nd);

%% Steady-state Kalman Gain

[~, L, ~] = kalman(kalman_sys,Qw,Rv);

%% initializing the system
%number of states, inputs and outputs
n_states  = length(sysd.A);
n_inputs  = size(sysd.B,2);
n_output  = size(sysd.C,1);

% Initial condition for the states and inputs
x1_initial= 0;
u1_initial= 0;
U         = [u1_initial];

%% calculating steady-state states and inputs
% Not necessary step for question but for validation of MPC solution

% Initial output set points are taken to be 0 and 0. 
y_sp1     = [y1_sp1];

% Steady-state targets for first step
ss_target1= [(eye(n_states)-A) -B;C zeros(n_output,n_inputs)]\[zeros(n_states,1);y_sp1];
x1_s      = ss_target1(1);
u1_s      = ss_target1(2);

% Unit step change for output set point is given after 25 mpc time steps
y_sp2     = [y1_sp2]';

%% setting up mpc

% State vector and input vector
x          = [x1_initial];
u          = [];

% Initial guess for the optimization variables
x1_0      = x1_initial*ones(1,N);
u1_0      = u1_initial*ones(1,N+1);
X0        = [x1_0 u1_0];

% Output definition
y         = C*x;
y_noise   = C*x + 0.1*(randn(1,1)-0.5);
E_plot    = [0];          % cummulative square error

%% MPC Simulation

U_delay = zeros(1,theta/del_t);

%mpc_steps = 5

for i       = 1:mpc_steps
    
    % Linear inequality constraint
    A_ineq  = [];
    b_ineq  = [];
    
    % Linear equality constraint
    Aeq     = [];
    beq     = [];
    
    % Bounds definition 
    lb_x1   = -inf*ones(1,N);   % Lower bound for state x1
    lb_u1   = -inf*ones(1,N+1); % Lower bound for control input u1
    
    ub_x1   = inf*ones(1,N);    % Upper bound for state x1
    ub_u1   = inf*ones(1,N+1);  % Upper bound for control variable u1

    
    % Lower and upper bounds for all the optimization variables
    lb      = [lb_x1 lb_u1];
    ub      = [ub_x1 ub_u1];
    
    % To change the steady state targets as per output set point after 25
    % time steps
    if i < step_change
        ysp = y_sp1;
    else
        ysp = y_sp2;
    end
    
    % fmincon optimization options
    options = optimoptions('fmincon','Display','none','Algorithm','interior-point','MaxFunctionEvaluations',4e3);
    nonlcon = @nonlinear_constraint;
    call_fun= @objective_function;
   
    % fmincon NLP solver
    % opt_var - optimized variable output, fval - function value
    [opt_var, fval] = fmincon(@(opt)call_fun(opt,C,N,ysp),...
    X0,A_ineq,b_ineq,Aeq,beq,lb,ub,@(opt)nonlcon(opt,A,B,x1_initial,U_delay,N),options);
    
    % Initial guess for optimization variables for the next prediction horizon
    X0      = opt_var;
    
    % take (time_delay*time_step +1)th input from opt
    U_1     = opt_var((1*N+theta/del_t+1));
    % add to previous input sequence
    U_delay = [U_delay U_1];
    
    % First input for system simulation
    U = U_delay(1);
    
    % get new initial first 5 U for next loop
    U_delay = U_delay(2:end);
    
    %Plant model simulation with first step input
    % Simulating the plant model for single time step with control input U0 from fmincon optimization
    mpc_sim = A_tilda*([x1_initial;U(1);0])....
        + B_tilda*(U) + sqrt(Qw)*(rand(na+nd,1)-ones(na+nd,1)*0.5);
    
    y_mpc   = C_tilda*mpc_sim + sqrt(Rv)*(rand(p,1)-ones(p,1)*0.5);
       
    %Estimating the state and disturbance using observer
    x_est   = A_tilda*([x1_initial;U(1);mpc_sim(3)])....
        + B_tilda*(U) + L*(y_mpc - C_tilda*([x1_initial;U(1);mpc_sim(3)]));
    
    %Steady-state target selection
    ss_l    = [eye(1)-A,-B;C,zeros(1)];
    ss_r    = [B_disturbance*x_est(3);ysp-Cd*x_est(3)];
    x_u_ss  = ss_l\ss_r;
    
    % Steady state targets
    x1_s    = x_u_ss (1);
    u1_s    = x_u_ss (2);

    %Initial state condition for next prediction horizon
    x1_initial = x_est(1);
    
    % Saving the output states for eeach time step
    x       = [x x_est(1)];
    
    % Saving the outputs for each time step
    y_noise = [y_noise y_mpc];
    y       = [y C*x_est(1)];
    E_plot  = [E_plot (E_plot(end) + (y(end) - ysp)^2)];
    
    % Concatanating the mpc control inputs
    u       = [u U_1];
end
toc

% %% Plotting the results
T         = t_0:del_t:mpc_steps*del_t;

% Plot for the state outputs
figure;plot(T,x,'b','Linewidth',2);grid on
xlabel('Simulation Time in seconds');ylabel('States')
title('System State Trajectory');legend('X1 - State')

% Plot for the control input
figure();stairs(T(1:end-1),u,'b','Linewidth',2);grid on
xlabel('Simulation Time in seconds');ylabel('Control Input')
title('System Control Input');legend('Control Input - u1')

% Plot for the outputs
y1_s    = [y1_sp1*ones(1,step_change) y1_sp2*ones(1,length(T)-step_change)];
y1_s(1) = 0;

figure();stairs(T,y1_s,'g--','Linewidth',2);grid on;hold on
plot(T,y(1,:),'b','Linewidth',2);hold on
plot(T,y_noise(1,:),'r','Linewidth',1)
xlabel('Simulation Time in seconds');ylabel('Outputs')
title('System Output Trajectory')
legend('Output Set Point - y1','Filtered Output y1','Unfiltered Output y1')

% Plot Integral square error
figure()
plot(T,E_plot,'r-','Linewidth',2);xlabel('Time [sec]') 
grid on ;legend('Integral square error'),title('Integral square error')
ylabel('Integral square error');

%% The end