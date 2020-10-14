function [opt_fun] = objective_function(opt,C,N,ysp)
C;
Q = eye(1);
R = eye(1);
S = eye(1);
opt_fun = 0;

% cost = cost + output regulation penalty + input rate regulation penalty 

for i = 1:N
    opt_fun = opt_fun + (([opt(i)]*C - [ysp])')*Q*([opt(i)]*C - [ysp])...
        + (([opt(i+1*N+1)] - [opt(i+1*N)]))*S*([opt(i+1*N+1)] - [opt(i+1*N)]);
    
    % input penalty (not included)
    %+ (([opt(i+1*N)] - [u1_s])')*R*([opt(i+1*N)] - [u1_s])...;
end

end