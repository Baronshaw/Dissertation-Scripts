function y = sinehole( ar1, Cr1, x )
        y = Cr1*sin(pi.*x./ar1)./(pi.*x/ar1);
        y( x==0 ) = Cr1;
end