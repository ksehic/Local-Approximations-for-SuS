 function f=force(t,x0)
 
S = 0.03; % [m^2/s^3] intensity of white noise

w = linspace(1*30*pi/300,15*pi,length(x0)/2);

delta_w = w(2)-w(1);

sigma_S = sqrt(2*S*delta_w);
 
f = -sigma_S*sum(x0(1:end/2).*cos(w.*t) + x0(end/2+1:end).*sin(w.*t),2);

 end