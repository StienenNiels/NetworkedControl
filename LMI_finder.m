function [Fprox,Gprox] = LMI_finder(h,alpha1,alpha2)
%LMI_finder
%    [Fprox,Gprox] = LMI_finder(H,ALPHA1,ALPHA2)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    02-Jun-2024 12:15:08

t2 = exp(h);
t3 = h.*(4.7e+1./1.0e+1);
t4 = -t3;
t6 = t2.*(8.5e+1./5.7e+1);
t5 = exp(t4);
t7 = -t6;
Fprox = reshape([t5,0.0,0.0,t5.*(8.5e+1./5.7e+1)+t7,t2,0.0,alpha1.*3.172825681224337e-1+alpha2.*(8.5e+1./5.7e+1)-t5.*3.172825681224337e-1+t7,-alpha2+t2,0.0],[3,3]);
if nargout > 1
    Gprox = [alpha2.*(-8.5e+1./4.7e+1)+8.5e+1./4.7e+1;alpha2-1.0;1.0];
end
end
