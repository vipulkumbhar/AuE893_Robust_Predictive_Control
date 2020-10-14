function [opt_fun] = obj_fun_c(opt,X1_set,X2_set,Usp1,Usp2,Q,R,Hp)

% Constatnts in the objective function
opt_fun_temp = 0;
U1           = Usp1;
U2           = Usp2;

% opt(i)       is x1
% opt(i+Hp)    is x2
% opt(i+2Hp)   is U1
% opt(i+1+3Hp) is U2

% Cost penalizing the tracking error
% Obj_function = Xhat'QXhat + Uhat'*R*Uhat + terminal cost

for i = 1:Hp
      
    cost_x           = ([opt(i)-X1_set opt(i+Hp)-X2_set])*Q*([opt(i)-X1_set;opt(i+Hp)-X2_set]);
    cost_u           = ([opt(i+2*Hp)-U1;opt(i+3*Hp+1)-U2])'*R*([opt(i+2*Hp)-U1;opt(i+3*Hp+1)-U2]);
    cost_delta_u     = 0;
    opt_fun_temp     = opt_fun_temp + 0.5*(cost_x +cost_u+cost_delta_u);
    
end

% Total cost function along with the terminal cost 
opt_fun = opt_fun_temp; %