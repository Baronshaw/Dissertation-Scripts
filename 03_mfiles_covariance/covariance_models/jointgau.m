function [z] = jointgau(ar1,at1,Cr1,x,y)
% this is the function for the joint gaussian covariance model

    z = Cr1.*( exp(-(sqrt(3)*x./ar1).^2).*exp(-(sqrt(3)*y./at1).^2) );

end