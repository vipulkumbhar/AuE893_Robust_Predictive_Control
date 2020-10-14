%% AuE 893: Robust Predictive Conrtol HW03, Excercise 1.61
% Author: Vipul Kumbhar
% ref:https://ijireeice.com/wp-content/uploads/2013/03/IJIREEICE2E_S_diwakar_Different_PID.pdf

clc;clear;close all;tic
t_end          = 80;
time_step      = 1;

s              = tf('s');
theta          = 5;
k              = 1;
tou            = 1;

% Ref:https://www.mathworks.com/help/mpc/ug/plant-models-with-delays.html
G              = tf([k],[tou 1],'IOdelay',theta,'TimeUnit','seconds');
tf_continuous  = G;
tf_descrete    = c2d(G,time_step)
sys            = ss(tf_descrete);

% State space model
A              = sys.A;
B              = sys.B;
C              = sys.C;
D              = sys.D;

% Check system response
u_sim          = 1*ones(1,(t_end/time_step +1));
t_sim          = 0:time_step:t_end;
x_ini          = [0];

% ODE45 equivalent, system response
[y,t,x]        = lsim(sys,u_sim,t_sim,x_ini);

% Plot system response for initial simulation
figure('Name','System Response');subplot(3,1,1);plot(t,y)
title('Output vs time plot')    ;legend('Y1')

subplot(3,1,2);plot(t,x);title('States vs time plot')
legend('X1')

subplot(3,1,3);plot(t,u_sim);title('Input vs time plot');legend('U1')
xlabel('Time (sec)')
close all
%% ZIEGLER-NICHOLS METHOD for PID tuning
% The Ziegler-Nichols design methods are the most popular methods used in
% process control to determine the parameters of a PID controller 

% The unit step response method is based on the open-loop step response of 
% the system.The unit step response of the process is characterized by two
% parameters, delay L1 and time constant T.
% These are determined by drawing a tangent line at the inflexion point, 
% where the slope of the step response has its maximum value. The 
% intersections of the tangent and the coordinate axes give the process 
% parameters.

[step_y,step_t,step_x] = step(sys);
slope                  = zeros(length(step_y),1);

for i = 1:length(step_y)-1
    slope(i,1)     = step_y(i+1) - step_y(i); 
end

[Max_slope, Index] = max(slope);
Max_y = step_y(Index+1);
Max_t = step_t(Index+1);

ZN_T  = 1/ Max_slope;                             % Delay
ZN_L1 = Max_t - (step_y(Index+1) - 0)/Max_slope;  % Time constant

% Kp, Ki, Kd parameters by ZIEGLER-NICHOLS method
Kp    = 1.2*ZN_T/ZN_L1
Ti    = 2*ZN_L1;
Td    = 0.1*ZN_L1;
Kd    = Td
Ki    = 1/Ti

%%
Y_setpoint     = 1;
U              = [0 0];
x_ini          = 0;
e_int          = 0;
e_prev         = 0;
E_plot         = [0];

Y_ref = [0];

for i=1:t_end/time_step
    t_sim      = 0:time_step:(time_step*i);

    % ODE45 equivalent, system response
    [y,t,x]    = lsim(sys,U,t_sim,x_ini);
    e          = Y_setpoint - y(end);
    
    % PID controller with ZN method for parameters
    u          = Kp*e + Ki*(e+e_int) + Kd*(e-e_prev); 
    
    % Futher optimized PID parameters
    u          = 0.45*e + 0.097*(e+e_int) +0.1*(e-e_prev);
    
    e_prev     = e;
    e_int      = e_int + e;
    U          = [U u];
    E_plot     = [E_plot (E_plot(end)+e*e)];
    
    % Check response for varied setpoints
    if i == 90
        Y_setpoint = 4;
    end
    if i== 180
        Y_setpoint = 0;
    end
    
    % Setpoints reference
    Y_ref = [Y_ref Y_setpoint];  
end
toc

figure()
plot(t_sim,y,'b','Linewidth',2); xlabel('Time [sec]')     
grid on ;title('System Output Trajectory')
hold on; ylabel('Outputs')
stairs(t_sim,Y_ref,'g--','Linewidth',2); 
legend('Output y1','Output Set Point - y1')

figure();stairs(t_sim,U(2:end),'b','Linewidth',2);
xlabel('Time [sec]'); ylabel('Control Input') 
grid on;legend('Control Input - u1');title('System Control Input')

figure()
plot(t_sim,x,'b','Linewidth',2);xlabel('Time [sec]'),ylabel('States')         
grid on ;legend('X1 - State'),title('System State Trajectory')

figure()
plot(t_sim,E_plot,'r-','Linewidth',2);xlabel('Time [sec]') 
grid on ;legend('Integral square error'),title('Integral square error')
ylabel('Integral square error');

%% close all