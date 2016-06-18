function [] = revisit_softdata()
% created 4/15/2014, this is creating figures that revisit the soft data to
% improve lambda1 and lambda2

deltaT = 365;
modplots = 0:5:50;
n = 3;
numbins = 10;
minpnts = 150; 
negval = 0;

% first, plot the CTM data, lambda1, and lambda2 every day for 2001
dayz = datevec(datenum(2001,1,1):datenum(2001,12,31));
dayz = datevec(datenum(2005,7,1):datenum(2005,12,31));
dayzdisp = dayz(:,1).*10000 + dayz(:,2).*100 + dayz(:,3);
[r c] = size(dayz);

for i = 1:r
    
    disp(dayz(i,1:3));
    plotCTM(dayz(i,1:3)); 
    plotLambdas(dayz(i,1:3),dayzdisp(i),n,deltaT,numbins,minpnts,modplots,negval);
    
    if length(findall(0,'type','figure')) >= 20, close all; end
    
end

end