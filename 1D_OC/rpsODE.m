% Phage resistant phage saturation model ODE with phage treatment control
function dxdt = rpsODE(t,x,p,u)
dxdt = [0;0;0;0]; % (S,R,P,I) -> (rescale) -> (x1,x2,x3,x4)

%----------- Loading parameters-----------
% susceptible bacteria growth rate
r = p.r;
% resistant bacteria growth rate
rp = p.rp;
% total bacteria carrying capacity
Kc = p.Kc;
% adsorption rate of phage:
phi = p.phi;
% phage density at half saturation
Pc = p.Pc;
% immune response killing rate parameter:
ep = p.ep;
% bacterial conc. at which immune response is half as effective:
Kd = p.Kd;
% burst size of phage:
beta = p.beta;
% decay rate of phage:
w = p.w;
% maximum growth rate of immune response:
a = p.a;
% max capacity of immune response:
Ki = p.Ki;
% conc. of bacteria at which imm resp growth rate is half its maximum:
Kn = p.Kn;
% probability of emergence of phage-resistant mutation per cell division
m = p.m;
% phage injection amount is accomplished in 1hr
P_in = p.P_in; 


%----------- rescaled parameters-----------
% recaled parameters will be noted by '(.)0' notation
ep0 = ep*Ki; 
q0 = P_in/Pc; 
Kcd0 = Kc/Kd; 
Knd0 = Kn/Kd; 
Kpd0 = Pc/Kd;
phi0 = phi*Kd;
% infection rate function (PS-R model)
F_P = phi0*(x(3)/(1 + x(3)));


n = x(1) + x(2);  % total (recaled) bacteria
dxdt(1) = r*x(1)*(1 - n/Kcd0)*(1-m) - Kpd0*x(1)*F_P - ep0*x(1)*x(4)/(1 + n);
dxdt(2) = rp*x(2)*(1 - n/Kcd0) + r*x(1)*(1 - n/Kcd0)*m - ep0*x(2)*x(4)/(1 + n);
dxdt(3) = (beta*x(1)*F_P) - phi0*x(1)*x(3) - w*x(3) + q0*u(1);
dxdt(4) = a*x(4)*(1 - x(4))*(n/(n + Knd0));
end

