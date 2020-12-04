function [opt_fun] = objective_function(opt,Y_setpoint,U_setpoint,N,Q,R)
opt_fun = 0;

for i = 1:N
    opt_fun = opt_fun + (([opt(i); opt(i+N)] - Y_setpoint)')*Q*([opt(i);opt(i+N)] - Y_setpoint)...
        + ((opt(i+2*N) - U_setpoint)')*R*(opt(i+2*N) - U_setpoint);
        
end
end