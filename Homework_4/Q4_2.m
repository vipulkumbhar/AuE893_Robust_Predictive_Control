%% AuE 893: Robust Predictive Conrtol HW04, Question 4
% Author: Vipul Kumbhar
% ref:

clc;clear;close all

% Given state space and gains
A     = [1 1;0 1];     B     = [0;1];
W     = [-0.1 0.1];    K     = [-0.4,-1.2];
Ak    = A + B*K;

% Equality constraints
% Cj*X <= Dj
Cj    = [1;0];
Dj    = [2;0];
% Aj*U <= Bj
Aj    = [1;-1];
Bj    = [1;1];

% Steps
N     = 50;
alpha = 0.0001;

% Disturbances
wmax  = 0.1;
W_a   = [wmax, wmax, -wmax, -wmax;
     wmax, -wmax, -wmax, wmax];

% Initialize parameters
Sum_AkW          = 0;
state_constraint = zeros(2,51);
input_constraint = zeros(2,51);

% Bounding tube by closed loop control sequence (Ak = A+BK)
for i =1:50+1
    Sum_AkW = Sum_AkW + (Ak^(i-1))*W_a;
    
    psi     = max(max(Cj(1)*Sum_AkW));   % X_u
    psi2    = max(max(Cj(2)*Sum_AkW));   % X_l
    
    theta   = max(max(Aj(1)*K*Sum_AkW)); % U_u
    theta2  = max(max(Aj(2)*K*Sum_AkW)); % U_u
    
    state_constraint(1,i) = (Dj(1)-(1-alpha)^-1*psi)/Cj(1);
    state_constraint(2,i) = (Dj(2)-(1-alpha)^-1*psi2)/Cj(2);
    input_constraint(1,i) = (Bj(1)-(1-alpha)^-1*theta)/Aj(1);
    input_constraint(2,i) = (Bj(2)-(1-alpha)^-1*theta2)/Aj(2);
end

N_lin   = 0:N;
X_u_lin = Dj(1)*ones(length(N_lin));
U_u_lin = Bj(1)*ones(length(N_lin));
U_l_lin = -Bj(2)*ones(length(N_lin));

figure('name','Final tightened constraints')
plot(N_lin,state_constraint(1,:),'r',...
    N_lin,input_constraint(1,:),'b',...
    N_lin,input_constraint(2,:),'g',...
    N_lin,X_u_lin,'r--',...
    N_lin,U_u_lin,'b--',...
    N_lin,U_l_lin,'g--',...
    'LineWidth',2)
legend('X_1 upper tightened constraint', 'U upper tightened constraint',...
    'U lower tightened constraint','X_1 Upper limit','U upper limit',...
    'U lower limit');grid on
%axis([1 50 -1.5 3.5])

%% Q4
pd = makedist('uniform')

X_startpoints_all = [1 1 0 -1;1 0 1 -1];
f_val_all =[];
W_all     = (rand(2,N-10)-0.5)/5;

% W_all_1 = ([rand(1,12) 0.2*ones(1,N-10-12)] - 0.5)/5;
% W_all_2 = ([rand(1,17) -0.3*ones(1,N-10-17)] - 0.5)/5;
% W_all =[W_all_1;W_all_2]; 

Y_all     = [];

for i_s = 1:4
    X_initial = X_startpoints_all(:,i_s);
    A         = A;
    C         = eye(2);
    D         = [0;0];
    Q         = 1*eye(2);
    R         = 1;

    Y_initial = X_initial;
    Y_setpoint= [0;0];             % X_setpoint

    U_setpoint= 0;  %(eye(2) - A)*Y_setpoint\B;

    time_step = 0.5;
    N_mpc_all = N;
    Tph       = 5;                 % [sec] prediction horizon
    t_mins    = Tph*2;
    N_mpc     = Tph/time_step;

    % Initial guesses
    U_guess = zeros(1,N_mpc+1);
    X_guess = X_initial.*ones(2,N_mpc);
    XU_guess= [X_guess(1,:) X_guess(2,:) U_guess];

    % Initialize outputs data for plots
    u     = 0;
    y     = X_initial;
    e     = [0;0];
    
    for i=1:N_mpc_all-t_mins
    
        % Linear inequality constraint
        A_ineq = [];
        b_ineq = [];
    
        % Linear equality constraint
        Aeq    = [];
        beq    = [];
        
        % Bounds definition
        lb_x1  = -inf*ones(1,N_mpc);     % Lower bound for state x1
        lb_x2  = -inf*ones(1,N_mpc);     % Lower bound for state x2
        lb_u1  = input_constraint(2,i:i+N_mpc) -  K*e;% Lower bound for control input u1
    
        ub_x1  = state_constraint(1,i:i+N_mpc-1);% Upper bound for state x1
        ub_x2  = inf*ones(1,N_mpc);      % Upper bound for state x2
        ub_u1  = input_constraint(1,i:i+N_mpc)- K*e;% Upper bound for control input u1
    
        % Lower and upper bounds for all the optimization variables
        lb = [lb_x1 lb_x2 lb_u1];
        ub = [ub_x1 ub_x2 ub_u1];
    
        % fmincon optimization options
        options = optimoptions('fmincon','Display','none','Algorithm','interior-point',....
        'MaxFunctionEvaluations',30000);
        nonlcon = @nonlin;
        call_fun= @objective_function;
    
        % fmincon NLP solver
        % opt_var - optimized variable output, fval - function value
        [opt_var, fval] = fmincon(@(opt)call_fun(opt,Y_setpoint,U_setpoint,N_mpc,Q,R),...
          XU_guess,A_ineq,b_ineq,Aeq,beq,lb,ub,@(opt)nonlcon(opt,A,B,X_initial,N_mpc),options);
      
        f_val_all(i_s,i) = fval;
        
         % Initial guess for optimization variables for the next prediction horizon
         XU_guess = opt_var;
         
         U  = opt_var((2*N_mpc)+1);
         
         % Control plus feedback
         U_bar   = U + K*e;
         
         %W  = wmax*[(rand - 0.5)/0.5;(rand - 0.5)/0.5];   % W E [-0.1 0.1]
         W = W_all(:,i);
         
         % Simulate systme
         mpc_sim = A*X_initial + B*U_bar + W;
         
         % Predicted error
         e = Ak*e + W;
         
         X_initial = mpc_sim;
    
         % Save the outputs for eeach time step
         y = [y mpc_sim];
    
        % Concatanate the mpc control inputs with initial inputs
         u = [u U_bar];
        
    end

    N_lin = 0:N_mpc_all-t_mins;
    x1_setpoint = Y_setpoint(1)*ones(1,length(N_lin));
    x2_setpoint = Y_setpoint(2)*ones(1,length(N_lin));
    X1_tightened_constraint = state_constraint(1,1:length(N_lin));

    U_setpoint = U_setpoint*ones(1,length(N_lin));
    U_u_tightened_constraint = [input_constraint(1,1) input_constraint(1,1:length(N_lin)-1)];
    U_l_tightened_constraint = [input_constraint(2,1) input_constraint(2,1:length(N_lin)-1)];
    U_u = 1*ones(1,length(N_lin));
    U_l = -1*ones(1,length(N_lin));

    if i_s ==4
        close all
    end
    
    figure('name','Part C - states plot')
    plot(N_lin,y(1,:),'r',...
        N_lin,x1_setpoint,'r--',...
        N_lin,y(2,:),'b',...
        N_lin,x2_setpoint,'b--',...
        N_lin,X1_tightened_constraint,'g--','LineWidth',2)

    title('System State Trajectory')
    xlabel('Simulation Time in seconds')
    ylabel('States');grid on
    legend('X1 - State','X1 - Set Point','X2 - State','X2 - Set Point','X1 upper tightened constraint')
    %axis([0 N_mpc_all-10 -2.5 2.5])

    figure('name','Part C - input plot')
    stairs(N_lin,u,'r','LineWidth',2)
    hold on
    plot(N_lin,U_u_tightened_constraint,'b--',...
      N_lin,U_l_tightened_constraint,'g--',...
      N_lin,U_u,'b',...
      N_lin,U_l,'g',...
      N_lin,U_setpoint,'--',...
     'LineWidth',2)
    grid on;title('System Control Input')
    xlabel('Simulation Time in seconds');ylabel('Control Input');grid on
    %axis([0 N_mpc_all-10 -1.5 1.5])
    legend('U - control input','U - upper tightened constraint',...
        'U - lower tightened constraint','U - upper limit','U - lower limit',...
         'U - setpoint')
    
    Y_all(2*i_s-1,:) = y(1,:);
    Y_all(2*i_s,:)   = y(2,:);
end

% Plot all
Ysize = size(Y_all);

figure()
for i=1:Ysize(1)/2
    plot(Y_all(2*(i-1)+1,:),Y_all(2*i,:),'LineWidth',2)
    hold on
end
title('System State Trajectory')
xlabel('X1')
ylabel('X2');grid on
axis([-3 3 -1.5 1.5])
lgd = legend('States [1 1]','States [1 0]','States [0 1]','States [-1 -1]');
title(lgd,'Initial state points = ')

% Plot disturbance
figure()
plot(W_all(1,:),'LineWidth',2);hold on
plot(W_all(2,:),'LineWidth',2);grid on
xlabel('N');ylabel('value')
legend('W1','W2')
title('State disturbance sequence')

%% The end