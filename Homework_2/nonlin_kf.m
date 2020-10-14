% Function to implement the nonlinear constraint to the optimization

function [c,ceq] = nonlin_kf(opt,X_ini,A,B,time_step,Hp)

% constants in the discretized nonlinear constraint function
A;
B;
time_step;
% Non-linear inequality constraint 
c = [];
I = eye(3);

ceq = [];

for i = 1:Hp
    
    % State constraint
    % 0         = x(k+1) - (I-time_step*A)*X(k) - time_step*B*U(k)
    
    % For descrete model
    ceq_temp    = [opt(i) opt(i+Hp) opt(i+2*Hp)]' - (A)*[X_ini(1) X_ini(2) X_ini(3)]' - B*[opt(i+3*Hp) opt(i+1+4*Hp)]';
    
    ceq(i)      = ceq_temp(1);
    ceq(i+Hp)   = ceq_temp(2);
    ceq(i+2*Hp) = ceq_temp(3);
    
    % For continuous model
    %ceq         = [opt(i) opt(i+Hp) opt(i+2*Hp)]' - (I-time_step*A)*[X_ini(1) X_ini(2) X_ini(3)]' - time_step*B*[opt(i+3*Hp) opt(i+1+4*Hp)]';
     
    % for next step, update xk to be x(k+1) of previous step
    X_ini       = [opt(i) opt(i+Hp) opt(i+2*Hp)]';
    
end