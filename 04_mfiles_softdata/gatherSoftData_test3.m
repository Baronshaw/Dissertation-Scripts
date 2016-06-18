function [] = gatherSoftData_test3()
% this function will load all the soft data calculated and extract only the
% variables that are needed for the BME analysis

n = 3; 
deltaT = 365; 
numbins = 10; 
minpnts = 150; 
negval = 0;

for j = [2001 2002]
    
    temp = datevec(datenum(j,1,1):datenum(j,12,31));
    dayWHI = temp(:,1:3);
    dayWHIdisps = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);
    [r c] = size(dayWHI);
    css = cell(r,1);
    lambda1 = cell(r,1);
    lambda2 = cell(r,1);
    limi = cell(r,1);
    probdens = cell(r,1);
    
    for i = 1:r
        disp([j i]);
        load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d_test3.mat',dayWHIdisps(i),deltaT,n,minpnts,numbins,negval));
        css{i,1} = [CTMlocs repmat(datenum(dayWHI(i,1),dayWHI(i,2),dayWHI(i,3)),length(CTMlocs),1)];
        lambda1{i,1} = meanGivMod(:,1);
        lambda2{i,1} = abs(varGivMod(:,1)); % no there's no zeros
        len = length(lambda1{i,1});
        if len > 0
            limi{i,1} = cell2mat(arrayfun(@linspace,lambda1{i,1}-1.96.*sqrt(lambda2{i,1}), ...
                lambda1{i,1}+1.96.*sqrt(lambda2{i,1}),13*ones(len,1),'uni',0));

            tempa =  mat2cell(limi{i,1},ones(size(limi{i,1},1),1),13);
            tempb = mat2cell(lambda1{i,1},ones(size(lambda1{i,1},1),1),1);
            tempc = mat2cell(lambda2{i,1},ones(size(lambda2{i,1},1),1),1);
            temp = cell2mat(cellfun(@truncNorm,tempa,tempb,tempc,'uni',0));
            probdens{i,1} = temp(:,2:end);  
        end
    end
    
    css = cell2mat(css);
    lambda1 = cell2mat(lambda1);
    lambda2 = cell2mat(lambda2);
    limi = cell2mat(limi);
    probdens = cell2mat(probdens);
    
    save(sprintf('../matfiles/PM2p5_soft_yr%d_test3.mat',j),'css','lambda1','lambda2','limi','probdens');

end

end