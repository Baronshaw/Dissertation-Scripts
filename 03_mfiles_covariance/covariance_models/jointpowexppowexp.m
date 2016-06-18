function [z] = jointpowexppowexp(alp,alpr1,alpt1,alpt2,alpr2,ar1,ar2,at1,at2,Cr1,x,y)
% this is the function for the joint double powered exponential covariance model

    z = Cr1*( (alp).*exp(-(3*x./ar1)).^alpr1.*exp(-(3*y./at1)).^alpt1 + ...
        (1-alp).*exp(-(3*x./ar2)).^alpr2.*exp(-(3*y./at2)).^alpt2 );

end