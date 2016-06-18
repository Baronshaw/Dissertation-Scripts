function [] = gatherMeanMod()
% this function will load all the mean modeled data calculated and extract only the
% variables that are needed for the RAMP simulation

n = 3; 
deltaT = 365; 
numbins = 10; 
minpnts = 150; 
negval = 0;

for j = 2001:2002
    
    temp = datevec(datenum(j,1,1):datenum(j,12,31));
    dayWHI = temp(:,1:3);
    dayWHIdisps = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);
    [r c] = size(dayWHI);
    css = cell(r,1);
    meanMod_all = cell(r,1);
    varMod_all = cell(r,1);
    limi = cell(r,1);
    probdens = cell(r,1);
    perctile_data_all = cell(r,1);
    
    for i = 1:r
        disp([j i]);
        load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat',dayWHIdisps(i),deltaT,n,minpnts,numbins,negval));
        css{i,1} = [cMSCTMv repmat(datenum(dayWHI(i,1),dayWHI(i,2),dayWHI(i,3)),length(dailyCTMg),1)];
        meanMod_all{i,1} = mean_Mod;
        varMod_all{i,1} = abs(var_Mod); % no there's no zeros
        perctile_data_all{i,1} = perctile_data;
    end
    
    css = cell2mat(css);
    meanMod_all = cell2mat(meanMod_all);
    varMod_all = cell2mat(varMod_all);
    perctile_data_all = cell2mat(perctile_data_all);

    save(sprintf('PM2p5_meanMod_yr%d.mat',j),'css','meanMod_all','varMod_all','perctile_data_all');

end

end