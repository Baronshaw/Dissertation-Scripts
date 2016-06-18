function [] = simulatelambdaGivMod2(yr2simu)
% this function will perform a simulation of RAMP and find lambda1 and
% lambda2 for given fixed modeled values from simulated lambda1s and
% lambda2s

if nargin < 1, yr2simu = 2001; end

% set seed
rng(0)

% get mean modeled for each percentile for each modeled space/time location
load(sprintf('PM2p5_meanMod_yr%d.mat',yr2simu));
% use css, meanMod_all

%%% generate lambda1 and lambda2 from a random field for each percentile of
%%% modeled values

% taken from: 
% http://www.mathworks.com/matlabcentral/fileexchange/27613-random-field-simulation
% -3/ar1 = -1/2c0 => c0 = ar/6

% get some covariance information
load('../matfiles/covmod_r_long_joint exponential exponential_joint.mat');
xdir = f.alp*f.ar1+(1-f.alp)*f.ar2;
ydir = f.alp*f.ar1+(1-f.alp)*f.ar2;
tdir = f.alp*f.at1+(1-f.alp)*f.at2;
corr.name = 'exp'; corr.sigma = f.Cr1; corr.c0 = [xdir/6 ydir/6 tdir/6]; corr.c1 = [];
idx = css(:,3) == datenum(2001,7,1);
cssSub = css(idx,:);
cd('gp/file_exchange3')
C = correlation_fun(corr,cssSub);
cd('../..')
corr.C = C; corr.A = []; corr.B = [];

% plotting each of the simulated lambda1 deciles for 2001/07/01
for i = 1:10
    disp(i);   
%     tic
%     cd('gp/file_exchange3')
%     [F{i},KL{i}] = randomfield(corr,cssSub);
%     cd('../..')
%     toc
    
    tic
    cd('gp/file_exchange3')
    [F{i},KL{i}] = randomfield(corr,cssSub,'filter',0.95);
    cd('../..')
    toc
    
    % save file
    save(sprintf('simulambda1_RF20010701_%0.2d.mat',i),'cssSub','corr','F','KL','-v7.3');
    
    % plot and save figure
    
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1
    
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(cssSub,F{i},lax,[plotax(:,1) plotax(:,2)],redpink);
    cax = [-14 14];
    cax = [-10 10];
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Randomly Generated Field of \\lambda_1 for Decile %0.2d',i));

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
        {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    for k = 1:length(allstates)
        cd ../09_mfiles_projections
        states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
        cd ../17_mfiles_validatelambda1
        plot(states(:,1),states(:,2),'k-');
    end 

    % save figure
    set(gcf,'Position',[0 0 800 500]);       
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 500]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/simulambda1_RF%0.2d.png',i));
      
end
close all

% generate lambda1 and lambda2 for all percentiles
lambda1_simu = normrnd(10,3,length(css),10);
lambda2_simu = normrnd(15,3,length(css),10);
% the s/t locations for these values are defined above in 'css'

%%%

% getting paired modeled data and space/time locations of obs
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yr2simu));
% use coordObs, yrmodaObs, Mod
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
dnyrmoda = datenum(yr,mo,da);
uniyrmoda = unique(dnyrmoda);

% calculate observed data where z_hat ~ N(lambda1,lambda2)
mObs_calc = NaN*ones(length(Mod),1);
vObs_calc = NaN*ones(length(Mod),1);
% 1) for every coordObs location, find closest css, find closest distCTMv
% 2) see Mod, see lambda1_simu, see lambda2_simu
% 3) calculate lambda1, calculate lambda2 through interpolation
for i = 1:length(uniyrmoda)
    disp(i);
    % subsetting all the gathered modeled data
    idx1 = uniyrmoda(i) == css(:,3); 
    cssSub = css(idx1,:); 
    meanMod_allSub = meanMod_all(idx1,:); 
    lambda1_simuSub = lambda1_simu(idx1,:);
    lambda2_simuSub = lambda2_simu(idx1,:);
    
    % subsetting paired data
    idx2 = uniyrmoda(i) == dnyrmoda;
    dnyrmodaSub = dnyrmoda(idx2);
    coordObsSub = coordObs(idx2,:);
    ModSub = Mod(idx2);
    
    % find closest points
    [idx, dist] = knnsearch(cssSub(:,1:2),coordObsSub);
    
    % calculate lambda1, lambda2 through interpolation
    % unfortunately all the inputs for 'interp1' can't be a matrix, so I
    % had to loop through each row
    mObs_interp = NaN*ones(length(idx),1); vObs_interp = NaN*ones(length(idx),1);
    x = meanMod_allSub(idx,:); y = lambda1_simuSub(idx,:); z = lambda2_simuSub(idx,:);
    for j = 1:length(ModSub)
        mObs_interp(j) = interp1(x(j,:),y(j,:),ModSub(j),'linear','extrap');
        vObs_interp(j) = interp1(x(j,:),z(j,:),ModSub(j),'linear','extrap');
    end
    %mObs_interp = interp1(meanMod_allSub(idx,:),lambda1_simuSub(idx,:),ModSub,'linear','extrap');
    %vObs_interp = interp1(meanMod_allSub(idx,:),lambda2_simuSub(idx,:),ModSub,'linear','extrap');
    
    % putting final values in the correct order
    [aidx bidx] = ismember([coordObsSub repmat(uniyrmoda(i),length(coordObsSub),1)],[coordObs dnyrmoda],'rows');
    mObs_calc(bidx) = mObs_interp;
    vObs_calc(bidx) = vObs_interp;
        
end

% save results
save(sprintf('simulatelambdaGivMod_%d.mat',yr2simu),'coordObs','yrmodaObs','Mod','mObs_calc','vObs_calc', ...
    'css','lambda1_simu','lambda2_simu','meanMod_all','perctile_data_all');

end