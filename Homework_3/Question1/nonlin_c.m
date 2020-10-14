function [c,ceq] = nonlin_c(opt,A,B,x1_initial,x2_initial,Hp)

% constants in the discretized nonlinear constraint function
% Non-linear inequality constraint 
c   = [];
ceq = [];

for i = 1:Hp
    
    % State constraint
    % 0         = x(k+1) - (A)*X(k) - B*U(k)
    % For descrete model
    ceq_temp    = [opt(i) opt(i+Hp)]' - (A)*[x1_initial x2_initial]' - B*[opt(i+2*Hp) opt(i+1+3*Hp)]';
    
    ceq(i)      = ceq_temp(1);
    ceq(i+Hp)   = ceq_temp(2);
  
    % for next step, update xk to be x(k+1) of previous step
    x1_initial  = opt(i); 
    x2_initial  = opt(i+Hp);
end

% Terminal constraints 
ceq(Hp)   = opt(Hp);
ceq(2*Hp) = opt(2*Hp);


