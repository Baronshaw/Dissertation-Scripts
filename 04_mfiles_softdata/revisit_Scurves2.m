function [] = revisit_Scurves2(soft,constant,gauss)
% created 4/16/2014 this function will explore S-curves with the worst
% errors

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% colors
colors = [ 'r' ; 'g' ; 'b' ];

% load cross-valdiation results
load(sprintf('../05_mfiles_crossvalidation/Xval_10fold%s%s%s_results_test1.mat', ...
    softstr,constr,gaussstr));

% variables to work with
unisID = unique(ckall(:,1:2),'rows');
for i = 1:length(unisID)
    idx = ckall(:,1) == unisID(i,1) & ckall(:,2) == unisID(i,2);
    obs(i) = mean(zhall(idx));
end

% look at all the absolute errors
idxloc = ckall(:,1) >= -2.5*10^6 & ckall(:,1) <= 2.5*10^6 & ckall(:,2) >= -2*10^6 & ckall(:,2) <= 1.5*10^6;
figure; hold on;
abserrall = abs(zkall-zhall);
hist(abserrall,100);
title('all absolute errors');

% poorest absolute errors
ub = prctile(abserrall,99.9);
len = sum(abserrall >= ub);
figure; hold on;
hist(abserrall(abserrall>=ub),100);
title('absolute errors in the upper precentile');
unisIDub = unique(ckall(abserrall>=ub&idxloc,1:2),'rows'); % stations with worst offenders

% obs v pred for the stations with worst abs errors with worst abs erros highlighted
for i = 1:10 % I shouldn't look at all of them, that's too many
    figure; hold on;
    idx = unisIDub(i,1) == ckall(:,1) & unisIDub(i,2) == ckall(:,2);
    plot(zkall(idx),zhall(idx),'b.');
    plot([0 max(zkall(idx))],[0 max(zkall(idx))],'k--'); % 1-1 line
    title(sprintf('mod v pred for loc %0.0f km %0.0f km',unisIDub(i,1)./1000,unisIDub(i,2)./1000));
    xlabel('pred');
    ylabel('obs');
    idx = unisIDub(i,1) == ckall(:,1) & unisIDub(i,2) == ckall(:,2) & abserrall >= ub;
    zksub = zkall(idx); zhsub = zhall(idx);
    for j = 1:sum(idx)
        plot(zksub(j),zhsub(j),'o','Color',colors(j));
    end
end

% should I show maps where these stations are???

% for the first 10 worst stations, for each worst point, show 3 closest S
% curves for each worst abs error
deltaT = 365; n = 3; minpnts = 150; numbins = 10; negval = 0;

% S curves
for i = 1:10 % showing 10 stations with worst offenders
    
    idx = unisIDub(i,1) == ckall(:,1) & unisIDub(i,2) == ckall(:,2) & abserrall >= ub;
    zhsub = zhall(idx);
    baddaysperstation = datevec(ckall(idx,3));
    baddaysdisp = baddaysperstation(:,1).*10000 + ...
        baddaysperstation(:,2).*100 + baddaysperstation(:,3);
    
    for j = 1:length(baddaysdisp)
        
        % loading soft data: soft data NOT from test1
        load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat', ...
            baddaysdisp(j),deltaT,n,minpnts,numbins,negval));

        distz = sqrt( (CTMlocs(:,1)-unisIDub(i,1)).^2 + (CTMlocs(:,2)-unisIDub(i,2)).^2 );
        [sorteddistz idxMod] = sort(distz);
        
        for k = 1:1
        
            figure; hold on;
        
            idxMod = distz == sorteddistz(k);
            plot(valMod(idxMod,:),valObs(idxMod,:),'b.');
            plot([get(gca,'XLim')],[get(gca,'XLim')],'k--');
            xlabel('Modeled Values of PM_{2.5} (\mug/m^3)');
            ylabel('Observed Values of PM_{2.5} (\mug/m^3)'); 
            title(sprintf('PM_{2.5} (\\mug/m^3) on %d, \\DeltaT=%d days, n=%d\nx=%0.0f km, y=%0.0f km',...        
                baddaysdisp(j),deltaT,n,CTMlocs(idxMod,1)./1000,CTMlocs(idxMod,2)./1000));
            
            % circle high abs error
            valModsub = valMod(idxMod,:); valObssub = valObs(idxMod,:);
            idx = (valObssub - zhsub(j)) == min(abs(valObssub - zhsub(j)));
            if sum(idx) == 1
                plot(valModsub(idx),valObssub(idx),'o','Color',colors(j));
            elseif sum(idx) ~= 1
                blah = 5;
            end

            % add bins
            for l = 1:numbins+1
                plot([perctile_data(idxMod,l) perctile_data(idxMod,l)],get(gca,'YLim'),'r-');
            end

            % calculate means in each bin
            plot(mean_Mod(idxMod,:),mean_Obs(idxMod,:),'ro','MarkerFaceColor','r');

            % linear interpolation in each bin
            xlim = get(gca,'XLim');
            ylim = get(gca,'YLim');
            interpx = linspace(xlim(1),xlim(2));
            interpy = interp1(mean_Mod(idxMod,:),mean_Obs(idxMod,:),interpx,'linear','extrap'); 
            plot(interpx,interpy,'r-');
            
            % add modeled value
            plot([dailyCTMg(idxMod) dailyCTMg(idxMod)],[ylim(1) meanGivMod(idxMod,1)], ...
                'c--','LineWidth',2);
            plot([xlim(1) dailyCTMg(idxMod)],[meanGivMod(idxMod,1) meanGivMod(idxMod,1)], ...
                'c--','LineWidth',2);
        
        end
    
    end

end

end