function [f,g] = maxlikelihood3withgradient(k)
global dt;

p = -k(4)*k(1)*exp(k(1)*dt) - k(5)*k(2)*exp(k(2)*dt) - (1-k(4)-k(5))*k(3)*exp(k(3)*dt); 
f=  -sum(  log(   p  ) );

if nargout > 1
g = [ sum(k(4)*(k(1)*dt+1).*exp(k(1)*dt)./p); sum( k(5)*(k(2)*dt+1).*exp(k(2)*dt)./p); sum((1-k(4)-k(5))*(k(3)*dt+1)*exp(k(3)*dt)./p); sum((k(1)*exp(k(1)*dt)  - k(3)*exp(k(3)*dt))./p); sum((k(2)*exp(k(2)*dt)  - k(3)*exp(k(3)*dt))./p)]; 
end

