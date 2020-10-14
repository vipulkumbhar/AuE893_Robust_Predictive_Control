%% AuE 893: Robust Predictive Conrtol HW03, Excercise 2.6 part b
% Author: Vipul Kumbhar
% ref:

clc;clear;close all;tic

% Important parameters
% Change alpha [1 0.1 0.08 0.05 0.01] and check system response
alpha          = 0.25; 

Q              = alpha*eye(2);
R              = eye(2);
N              = 3;              % N
del_t          = 1;
t_0            = 0;              % MPC start time
t_sim          = 80;             % Total simulation time
Tp             = del_t*N;        % Prediction horizon length
mpc_steps      = (t_sim - t_0)/del_t; %number of steps
step_change    = 25;             % y1 set point step change after 25 steps

% Given descrete state space model 
A              = [2 1;0 2];
B              = eye(2);
C              = eye(2);
D              = zeros(2,2);

% Compute terminal weight and LOR control law
% [K,P]          = dlqr(A,B,Q,R);
% Xmax           = [inf;5];
% Xmin           = [-inf;-5];
% Umax           = [1;1];
% Umin           = [-1;-1];

%% parameter definition
% initial states and input definition
x1_initial    = 0;
x2_initial    = 0;
u1_initial    = 0;
u2_initial    = 0;

n_states = length(A);      % number of states
n_inputs = size(B,2);      % number of inputs
n_output = size(C,1);      % number of outputs

% calculating steady-state states and inputs
% Initial output set points are taken to be 0 and 0. 
y1_sp1 = 0;
y2_sp1 = 0;
y_sp1 = [y1_sp1 y2_sp1]';

% Steady-state targets for first 25 step
ss_target1 = [(eye(n_states)-A) -B;C zeros(n_output,n_inputs)]\[zeros(n_states,1);y_sp1];
x1_s1 = ss_target1(1);
x2_s1 = ss_target1(2);
u1_s1 = ss_target1(3);
u2_s1 = ss_target1(4);

% Unit step change for output set point is given after 25 mpc time steps
y1_sp2 = 1;
y2_sp2 = 0;
y_sp2 = [y1_sp2 y2_sp2]';

% Steady-state targets for rest of the simulation
ss_target2 = [(eye(n_states)-A) -B;C zeros(n_output,n_inputs)]\[zeros(n_states,1);y_sp2];
x1_s2 = ss_target2(1);
x2_s2 = ss_target2(2);
u1_s2 = ss_target2(3);
u2_s2 = ss_target2(4);

%% MPC Formulation
% State vector and input vector
x           = [x1_initial; x2_initial];
u           = [];

% Initial guess for the optimization variables
x1_0        = x1_initial*ones(1,N);
x2_0        = x2_initial*ones(1,N);
u1_0        = u1_initial*ones(1,N+1);
u2_0        = u2_initial*ones(1,N+1);
X0          = [x1_0 x2_0 u1_0 u2_0];

% Output definition
y           = C*x;

%% Running the mpc
for i = 1:mpc_steps
    
    % Linear inequality constraint
    A_ineq = [];
    b_ineq = [];
    
    % Linear equality constraint
    Aeq = [];
    beq = [];
    
    % Bounds definition
    lb_x1 = -inf*ones(1,N);     % Lower bound for state x1
    lb_x2 = -inf*ones(1,N);     % Lower bound for state x2
    lb_u1 = -1*ones(1,N+1);   % Lower bound for control input u1
    lb_u2 = -1*ones(1,N+1);   % Lower bound for control input u2
    
    ub_x1 = 5*ones(1,N);      % Upper bound for state x1
    ub_x2 = 5*ones(1,N);      % Upper bound for state x2
    ub_u1 = 1*ones(1,N+1);    % Upper bound for control input u1
    ub_u2 = 1*ones(1,N+1);    % Upper bound for control input u2
    
    % Lower and upper bounds for all the optimization variables
    lb = [lb_x1 lb_x2 lb_u1 lb_u2];
    ub = [ub_x1 ub_x2 ub_u1 ub_u2];
    
    % To change the steady state targets as per output set point after 25
    % time steps
    if i < step_change
        x1_s = x1_s1;x2_s = x2_s1;u1_s = u1_s1;u2_s = u2_s1;
    else
        x1_s = x1_s2;x2_s = x2_s2;u1_s = u1_s2;u2_s = u2_s2;
    end
        
    % fmincon optimization options
    options = optimoptions('fmincon','Display','none','Algorithm','interior-point',....
        'MaxFunctionEvaluations',30000);
    nonlcon  = @nonlin_a;
    call_fun = @obj_fun_a;
    
    % fmincon NLP solver
    % opt_var - optimized variable output, fval - function value
    [opt_var, fval] = fmincon(@(opt)call_fun(opt,x1_s,x2_s,u1_s,u2_s,Q,R,N),...
        X0,A_ineq,b_ineq,Aeq,beq,lb,ub,@(opt)nonlcon(opt,A,B,x1_initial,x2_initial,N),options);
    
    % Initial guess for optimization variables for the next prediction horizon
    X0 = opt_var;
    
    % First time step control input from the fmincon optimization
    U = [opt_var((2*N)+1); opt_var((3*N)+2)];
    
    % Simulating the plant model for single time step with control input U0 from fmincon optimization
    mpc_sim = A*([x1_initial;x2_initial]) + B*(U);
    
    % Initial state condition for next prediction horizon
    x1_initial = mpc_sim(1);
    x2_initial = mpc_sim(2);
    
    % Save the output states for each time step
    x = [x mpc_sim];
    
    % Save the outputs for eeach time step
    y = [y C*mpc_sim];
    
    % Concatanate the mpc control inputs with initial inputs
    u = [u U];
end
toc
%% Plotting the results

T = t_0:del_t:t_sim;
% steady state target vector for plotting purpose
x1_s = [x1_s1*ones(1,step_change) x1_s2*ones(1,length(T)-step_change)];
x2_s = [x2_s1*ones(1,step_change) x2_s2*ones(1,length(T)-step_change)];

% Plot for the state outputs
figure
plot(T,x(1,:),'b',T,x1_s,'--b',T,x(2,:),'r',T,x2_s,'--r','Linewidth',2)
xlabel('Simulation Time in seconds')
ylabel('States')
title('System State Trajectory')
legend({'X1 - State','X1 - Set Point','X2 - State','X2 - Set Point'},'Location','southwest')
grid on

% Plot for the control input
figure
stairs(T(1:end-1),u(1,:),'b','Linewidth',2);hold on
stairs(T(1:end-1),u(2,:),'r','Linewidth',2)
xlabel('Simulation Time in seconds');ylabel('Control Input')
title('System Control Input')
legend('Control Input - u1','Control Input - u2');grid on

% Plot for the outputs
y1_s = [y1_sp1*ones(1,step_change) y1_sp2*ones(1,length(T)-step_change)];
y2_s = [y2_sp1*ones(1,step_change) y2_sp2*ones(1,length(T)-step_change)];
figure
plot(T,y(1,:),'b',T,y1_s,'--b',T,y(2,:),'r',T,y2_s,'--r','Linewidth',2)
xlabel('Simulation Time in seconds')
ylabel('Outputs')
title('System Output Trajectory')
legend({'Output 1','Output Set Point - 1', 'Output 2','Output Set Point - 2'},'Location','northwest')
grid on