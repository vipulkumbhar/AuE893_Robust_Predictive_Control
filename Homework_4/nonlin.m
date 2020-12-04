function [c,ceq] = nonlin(opt,A,B,X_ini,N)

c = [];
x1_ini = X_ini(1);
x2_ini = X_ini(2);

for i = 1:N

    ceq(i)   = opt(i) - A(1,:)*([x1_ini,x2_ini]') - B(1,:)*(opt(2*N+i));
    ceq(i+N) = opt(i+N) - A(2,:)*([x1_ini,x2_ini]') - B(2,:)*(opt(2*N+i));
    
    x1_ini = opt(i);
    x2_ini = opt(i+N);

end
end