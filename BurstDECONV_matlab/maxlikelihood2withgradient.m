function [f,g] = maxlikelihood2withgradient(k)
global dt;

p = -k(3)*k(1)*exp(k(1)*dt) - (1-k(3))*k(2)*exp(k(2)*dt); 
f=  -sum(  log(   p  ) );

if nargout > 1
g = [ sum(k(3)*(k(1)*dt+1).*exp(k(1)*dt)./p); sum((1-k(3))*(1+k(2)*dt).*exp(k(2)*dt)./p); sum((k(1)*exp(k(1)*dt)  - k(2)*exp(k(2)*dt))./p)]; 
end

