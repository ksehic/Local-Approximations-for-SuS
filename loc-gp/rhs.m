function xdot=rhs(t,x,input)
    
    alpha = 0.1;

    uy = 0.04;

    m =  6*1e4; %[kg]mass

    k = 5*1e6; %[N/m] stifness

    xi = 0.05; % viscous damping ratio

    c = 2*m*xi*sqrt(k/m); % damping
    
    A = 1;
    
    n = 3;
    
    gamma=0.5;
    
    beta=0.5;

    xdot_1 = x(2);
    
    xdot_2 = -(c/m)*x(2) - (k/m)*(alpha*x(1)+(1-alpha)*uy*x(3)) + force(t,input);

    xdot_3 = A*x(2) - beta * abs(x(2)) * ( abs(x(3)) )^(n-1) * x(3) - gamma * x(2)*(abs(x(3)))^n;
    
    xdot = [xdot_1 ; xdot_2 ; xdot_3];

 end