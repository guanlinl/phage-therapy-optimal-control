% Input : the unconstrained minimizer of QP
% Output: contrained solution by projection operator
function ustar = proj_constrain_set(uhat)
% define some directional vectors
e1 = [1;0]; e2 = [0;1]; v1 = ones(2,1); vg = [-1;1];
xi = (dot(v1, uhat) - 1)/2; 
u = uhat;
% active sets 
   if dot(e1,u) >= 0 && dot(e2,u) >= 0 && dot(v1,u) - 1 <=0
      ustar = u; % uhat is in the constrained set
   elseif dot(v1,u) - 1 >= 0 && dot(vg,u) + 1 >= 0 && dot(vg,u) - 1 <=0
      ustar = u - xi.*v1; 
   elseif dot(e2,u) - 1 >=0 && dot(vg,u) - 1 >=0
       ustar = [0;1];
   elseif dot(e1,u) <=0 && dot(e2,u) >=0 && dot(e2,u) - 1 <=0
       ustar = [0;u(2)];
   elseif dot(e1,u) <=0 && dot(e2,u) <=0
       ustar = [0;0];
   elseif dot(e1,u) >=0 && dot(e1,u) - 1 <=0 && dot(e2,u) <=0
       ustar = [u(1);0];
   else
      ustar = [1;0];
   end
end




