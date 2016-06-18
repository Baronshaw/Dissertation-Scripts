function y = pentaspherical( ar1, Cr1, x )
        y = Cr1.*(1-((15/8).*(x./ar1).^1-(5/4).*(x./ar1).^3) ...
            +(3/8).*(x./ar1).^5);
        y( x>ar1 ) = 0;
end