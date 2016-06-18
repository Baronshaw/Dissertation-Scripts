function [] = lambdacompare()
% this function will create estimates comparing mean bias/standard error
% with the approximation of lambda1 and lambda2

% load data
load('matfiles/traditional_performance_grid.mat');  
valplot = { num ; mObs ; mMod ; mBias ; nBias ; nmBias ; fBias ; ...
    mErr ; nErr ; nmErr ; fErr ; R ; R2 ; sBias ; msBias ; ...
    rmsBias ; nrmsBias ; mDsBias ; m2DmsBias ; s2DmsBias ; beta1 ; vObs ; vMod };  
for i = 1:length(valplot)   
    cax{i} = [prctile(valplot{i}(:),5) prctile(valplot{i}(:),95)]; 
end
yr = floor(yrNday./10^4); uniyr = unique(yr);
mo = floor((yrNday - yr*10^4)./10^2);
da = yrNday - yr*10^4 - mo.*10^2;
yrznum = datevec(datenum(uniyr,[1 6],1));
[a b] = size(yrznum);
datetoload = yrznum(1,1).*10^4 + yrznum(1,2).*10^2 + yrznum(1,3);
load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',datetoload(1)));
c = length(0:5:50);

% load the overall S-curve
load('matfiles/globalSCurve.mat');
valplot = { numGivMod ; mObsGivMod ; mModGivMod ; mBiasGivMod ; nBiasGivMod ; nmBiasGivMod ; fBiasGivMod ; ...
    mErrGivMod ; nErrGivMod ; nmErrGivMod ; fErrGivMod ; RGivMod ; R2GivMod ; sBiasGivMod ; msBiasGivMod ; ...
    rmsBiasGivMod ; nrmsBiasGivMod ; mDsBiasGivMod ; m2DmsBiasGivMod ; s2DmsBiasGivMod ; beta1GivMod ; vObsGivMod ; vModGivMod }; 
strnm = { 'number of paired modeled and obs'; 'mean obs value' ; 'mean modeled value' ; ...
    'mean bias' ; 'normalized bias' ; 'normalized mean bias' ; 'fractional bias' ; 'mean error' ; ...
    'normalized error' ; 'normalized mean error' ; 'fractional error' ; 'correlation' ; 'correlation squared' ; ...
    'standard bias' ; 'mean squared bias' ; 'root mean squared bias' ; 'normalized root mean squared bias' ; ...
    'mean bias DIV standard bias' ; 'mean bias squared DIV mean squared bias' ; ...
    'variance of bias DIV mean squared bias' ; 'beta1' ; 'variance of obs' ; 'variance of mod' };

% create plots and gifs for the overall S-curve
for l = 1:length(strnm)
    for i = 1:a
        disp(i);
        idx = yrznum(i,1) == yrmodaCTMv(:,1) & yrznum(i,2) == yrmodaCTMv(:,2) & yrznum(i,3) == yrmodaCTMv(:,3);

        for j = 1:c

            figure; hold on;

            load('../09_mfiles_projections/USAcontiguous.mat');
            cd ../09_mfiles_projections
            plotax = ell2lambertcc([x,y],'whiproj2001');
            cd ../13_mfiles_modelperformance

            % plotting 
            lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
            [xg yg Zg] = plotField(distCTMv(idx,:),repmat(valplot{l}(j),sum(idx),1),lax,[plotax(:,1) plotax(:,2)]);
            caxis(cax{l});   
            colorbar;
            axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

            % setting axis
            set(gca,'XTickLabel',get(gca,'XTick')/1000);
            set(gca,'YTickLabel',get(gca,'YTick')/1000);
            xlabel('km');
            ylabel('km');

            % overlaying the states
            load('../09_mfiles_projections/USAstates5.mat');
            for k = 1:length(X)
                cd ../09_mfiles_projections
                states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
                cd ../13_mfiles_modelperformance
                plot(states(:,1),states(:,2),'k-');
            end

            % title 
            title(sprintf('Global: %s = %d of PM_{2.5} (\\mug/m^3) on %s',strnm{l},modplots(j),datestr(yrznum(i,:))));    

            % save figure 
            set(gcf,'Position',[0 0 800 600]);
            set(gcf,'PaperUnits','inches');    
            set(gcf,'PaperPosition',[0 0 800 600]./100);
            set(gcf,'PaperPositionMode','manual');
            print(gcf,'-painters','-dpng','-r600',sprintf('figures/Global%s%d_%s_grid.png',strnm{l},modplots(j),datestr(yrznum(i,:))));

            drawnow;
            frameI = getframe(gcf);
            imI = frame2im(frameI);
            [imindI,cmI] = rgb2ind(imI,256);
            outfileI = sprintf('figures/Global%s_%s_grid.gif',strnm{l},datestr(yrznum(i,:)));
            % On the first loop, create the file. In subsequent loops, append.
            if j == 1
                imwrite(imindI,cmI,outfileI,'gif','DelayTime',1,'loopcount',inf);
            else
                imwrite(imindI,cmI,outfileI,'gif','DelayTime',1,'writemode','append');
            end

        end
        close all;
    end
end

end