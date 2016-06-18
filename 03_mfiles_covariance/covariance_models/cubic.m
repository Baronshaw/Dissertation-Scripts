function y = cubic( ar1, Cr1, x )
        y = Cr1.*(1-(7.*(x./ar1).^2-(8.75).*(x./ar1).^3 ...
            -(3.5).*(x./ar1).^5)-(0.75).*(x./ar1).^7);
        y( x>ar1 ) = 0;
end