function [z] = jointpowexp( alpr1, alpt1, ar1, at1, Cr1, x, y )
% this is the function for the joint powered exponential covariance model

    z = Cr1.*( exp(-(3*x./ar1).^alpr1).*exp(-(3*y./at1).^alpt1) );

end