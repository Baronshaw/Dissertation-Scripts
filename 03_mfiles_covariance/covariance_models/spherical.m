function y = spherical( ar1, Cr1, x )
        y = Cr1.*(1-((3/2).*(x./ar1)-(1/2).*(x./ar1).^3));
        y( x>ar1 ) = 0;
end