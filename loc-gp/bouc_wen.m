function x = bouc_wen(x0)
x0=x0.';

t_start = 0;
t_end = 8; % final time in seconds.

dt=110;

time_span = linspace(t_start,t_end,dt);

u0 = [0 0 0];

[t ,x ]= ode45 ( @(t,x) rhs(t,x,x0) , time_span , u0 );

x=x(end,1);

return