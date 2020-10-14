function [opt_fun] = obj_fun_kf(opt,time_step,X1_set,X2_set,X3_set,Usp,Ysp,C,Hp)

% Constatnts in the objective function

Q = C'*C;    % C'C
S = eye(2);
R = eye(2);
time_step;

opt_fun_temp = 0;
U1           = Usp(1);
U2           = Usp(2);

% opt(i)     is x1
% opt(i+Hp)  is x2
% opt(i+2Hp) is x3
% opt(i+3Hp) is U1
% opt(i+1+4Hp) is U2

% Cost penalizing the tracking error
% Obj_function = Xhat'QXhat + Uhat'*R*Uhat + delta_U'*S*deltau 
% here xhat = (C*Xs + Cd*d - Ysp)

for i = 1:Hp
       
    cost_x           = ([opt(i)-X1_set opt(i+Hp)-X2_set opt(i+2*Hp)-X3_set])*Q*([opt(i)-X1_set;opt(i+Hp)-X2_set;opt(i+2*Hp)-X3_set]);
    cost_u           = ([opt(i+3*Hp)-U1;opt(i+4*Hp+1)-U2])'*R*([opt(i+3*Hp)-U1;opt(i+4*Hp+1)-U2]);
    
    opt_fun_temp     = opt_fun_temp+0.5*(cost_x+cost_u);
    
end

% Total cost function along with the terminal cost 
opt_fun = opt_fun_temp + 4*opt(5*Hp+3)^2;

