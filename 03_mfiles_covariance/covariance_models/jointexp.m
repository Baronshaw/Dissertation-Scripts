function [z] = jointexp( ar1, at1, Cr1, x, y )
% this is the function for the joint exponential covariance model

    z = Cr1.*( exp(-(3*x./ar1)).*exp(-(3*y./at1)) );

end