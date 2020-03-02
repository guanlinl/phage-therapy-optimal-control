function dfdx = Jacobian_mat(x1,x2,x3,x4,x5,p)
dfdx = zeros(5,5);
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
F_P1 = phi0*(x3/(1 + x3));
F_P2 = phi0*(x5/(1 + x5));


% construct Jacobian matrix dfdx, size 5 by 5
% row 1
dfdx(1,1) = ((1 - m)*r*(Kcd0 - 2*x1 - x2)/Kcd0) - ...
    (Kpd0*F_P1) - (ep0*(1 + x2)*x4/(1 + x1 + x2)^2);
dfdx(1,2) = ((m - 1)*r*x1/Kcd0) + (ep0*x1*x4/(1 + x1 + x2)^2);
dfdx(1,3) = -Kpd0*phi0*x1/((1 + x3)^2);
dfdx(1,4) = -ep0*x1/(1 + x1 + x2);

% row 2
dfdx(2,1) = (m*r*(Kcd0 - 2*x1 - x2)/Kcd0) - rp*x2/Kcd0 + ...
    (ep0*x2*x4/(1 + x1 + x2)^2);
dfdx(2,2) = -m*r*x1/Kcd0 + rp*(Kcd0 - x1 - 2*x2)/Kcd0 - ...
    ep0*(1 + x1)*x4/((1 + x1 + x2)^2) - (Kpd0*F_P2);
dfdx(2,4) = -ep0*x2/(1 + x1 + x2);
dfdx(2,5) = -Kpd0*phi0*x2/((1 + x5)^2); 

% row 3
dfdx(3,1) = beta*F_P1 - phi0*x3;
dfdx(3,3) = beta*phi0*x1/((1 + x3)^2) - w - phi0*x1;

% row 4
dfdx(4,1) = -a*Knd0*x4*(x4 - 1)/((x1 + x2 + Knd0)^2);
dfdx(4,2) = -a*Knd0*x4*(x4 - 1)/((x1 + x2 + Knd0)^2);
dfdx(4,4) = a*(1 - 2*x4)*(x1 + x2)/(x1 + x2 + Knd0);

% row 5
dfdx(5,2) = beta*F_P2 - phi0*x5; 
dfdx(5,5) = beta*phi0*x2/((1 + x5)^2) - w - phi0*x2;
    
end