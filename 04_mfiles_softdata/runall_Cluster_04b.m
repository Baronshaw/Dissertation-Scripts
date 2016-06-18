function [] = runall_Cluster_04b(SOFTIDX)
% this function will run all the data in the 04_mfiles_softdata folder

if nargin < 1, SOFTIDX = 8; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../04_mfiles_softdata

load('../matfiles/prepSoft.mat');

% variables
deltaT = 365;
modplots = 0:5:50;
uniyrsOff = [ 1998 1999 ; 1999 1999 ; 1999 2000 ; 2000 2000 ; 2000 2001 ; ...
    2001 2001 ; 2001 2002 ; 2002 2002 ; 2002 2003 ; 2003 2003 ; 2003 2004 ; ...
    2004 2004 ; 2004 2005 ; 2005 2005 ; 2005 2006 ; 2006 2006 ; 2006 2007 ; ...
    2007 2007 ; 2007 2008 ; 2008 2008 ; 2008 2009 ; 2009 2009 ; 2009 2010 ; ...
    2010 2010 ; 2010 2011 ; 2011 2011 ; 2011 2012 ];
n = 3;
numbins = 10;
minpnts = 150; 
negval = 0;

% % added later
% n = 6;
% deltaT = 180;
%
% % added later
% n = 12;
% deltaT = 90;
%
% % added later
% n = 6;
% deltaT = 90;
% numbins = 5;
%
% added later
% n = 1,2,3,4,5,6,7,8,9,10,11;deltaT=365,180,90,45
deltaT = 90;
n = 3;
minpnts = 50; 

% finding soft data all the days within this range
temp = datevec([datenum(2001,1,1):datenum(2002,12,31) datenum(2005,1,1):datenum(2007,12,31)]);
temp = datevec([datenum(2001,1,1):datenum(2002,12,31)]);
temp = datevec([datenum(2005,1,1):datenum(2005,12,31)]);
temp = datevec([datenum(2001,1,1):datenum(2001,12,31)]);
dayWHI = temp(:,1:3); 

lendiv = 10; % how many jobs there will be

dayWHIall = temp(:,1:3);
[r c] = size(dayWHIall);
len = 1:floor(r/lendiv):r;
if SOFTIDX == lendiv
    dayWHI = dayWHIall(len(SOFTIDX):r,:);
else
    dayWHI = dayWHIall(len(SOFTIDX):len(SOFTIDX+1)-1,:);   
end

dayWHIdisp = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);
[r c] = size(dayWHI);
minT = datevec(datenum(dayWHI)-floor(deltaT/2));
maxT = datevec(datenum(dayWHI)+floor(deltaT/2));

for i = 1:r
    modplotss = modplots;
    disp(dayWHI(i,:));
      
    % loading the correct years 
    idx = minT(i,1)==uniyrsOff(:,1) & maxT(i,1)==uniyrsOff(:,2);
    tic  
    
    [meanGivMod,varGivMod,valMod,valObs,perctile_data,numpnts,idxtME, ...
        idxcMS,modplotss,dailyCTMg,mean_Mod,mean_Obs,var_Mod,var_Obs,CTMlocs,cMSCTMv, ...
        cMSObsx,cMSObsy,tMEObs] ...
        = getSoftData(dayWHI(i,:),n,deltaT,numbins,minpnts,modplotss,negval,D{idx},Obsg{idx},...
        CTMg{idx},cMSObs{idx},tMEO{idx},cMSCTM{idx},tMECTM{idx},...
        distCTMv{idx},yrmodaCTMv{idx},dailyCTMv{idx});
    
%     [meanGivMod,varGivMod,valMod,valObs,numpnts,idxtME, ...
%         idxcMS,modplotss,dailyCTMg,CTMlocs,cMSCTMv, ...
%         cMSObsx,cMSObsy,tMEObs] ...
%         = getSoftData_LinReg(dayWHI(i,:),n,deltaT,numbins,minpnts,modplotss,negval,D{idx},Obsg{idx},...
%         CTMg{idx},cMSObs{idx},tMEO{idx},cMSCTM{idx},tMECTM{idx},...
%         distCTMv{idx},yrmodaCTMv{idx},dailyCTMv{idx});
    
    toc
    
%     save(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d_LinReg.mat',dayWHIdisp(i),deltaT,n,minpnts,numbins,negval),...
%         'meanGivMod','varGivMod','valMod','valObs','numpnts',...
%         'idxtME','idxcMS','modplots','dailyCTMg','dayWHIdisp','CTMlocs','cMSCTMv', ...
%         'cMSObsx','cMSObsy','tMEObs');
    
    
    save(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat',dayWHIdisp(i),deltaT,n,minpnts,numbins,negval),...
        'meanGivMod','varGivMod','valMod','valObs','perctile_data','numpnts',...
        'idxtME','idxcMS','modplots','dailyCTMg','dayWHIdisp','mean_Mod','mean_Obs','var_Mod','var_Obs','CTMlocs','cMSCTMv', ...
        'cMSObsx','cMSObsy','tMEObs');

end

% % visualizing soft data with optimize parameters
% temp = unique(dayWHI(:,1:2),'rows');
% plotmonths = [temp ones(size(temp,1),1)];
% plotmonthsdisp = plotmonths(:,1)*10^4 + plotmonths(:,2)*10^2 + plotmonths(:,3);
% for i = 28:length(plotmonths)
%     disp(plotmonthsdisp(i,:));
%     if i<28 | i>36
%         plotCTM(plotmonths(i,:)); 
%         plotLambdas(plotmonths(i,:),plotmonthsdisp(i,:),n,deltaT,numbins,minpnts,modplots,negval);
%         plotLambdas_givMod(plotmonths(i,:),plotmonthsdisp(i,:),n,deltaT,numbins,minpnts,modplots,negval);
%         plotScurve(plotmonths(i,:),plotmonthsdisp(i,:),n,deltaT,numbins,minpnts,modplots,negval);
%     end
%     close all
% end

end