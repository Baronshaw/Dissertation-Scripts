function [z] = jointexppowexp(alp,ar1,alpr2,ar2,at1,alpt2,at2,Cr1,x,y)
% this is the function for the joint exponential powered exponential covariance model

    z = Cr1*( (alp).*exp(-(3*x./ar1)).*exp(-(3*y./at1)) + ...
        (1-alp).*exp(-(3*x./ar2)).^alpr2.*exp(-(3*y./at2)).^alpt2 );

end