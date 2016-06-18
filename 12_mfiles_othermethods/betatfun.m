function [betat] = betatfun(betatm1,rho,eta)
    betat = rho*betatm1 + eta;
end