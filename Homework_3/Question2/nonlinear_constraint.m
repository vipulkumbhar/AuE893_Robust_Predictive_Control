function [c,ceq] = nonlinear_constraint(opt,A,B,x1_ini,U_delay,N)
c      = [];
ceq    = [];
U_delay; % Last 5 inputs

for i = 1:N
    % Equality constraint
    ceq(i)       = opt(i) - A(1,:)*([x1_ini]) - B(1,:)*([opt(1*N+i)]);
    x1_ini       = opt(i);
    
    % previous delay inputs as first new inputs
    if i<length(U_delay)+1
        ceq(1*N+i) = opt(1*N+i)-U_delay(i);
    end

end
