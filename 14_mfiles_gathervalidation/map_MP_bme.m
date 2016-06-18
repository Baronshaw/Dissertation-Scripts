function [] = map_MP_bme(forceddist,stationstr)
% create maps at the monitoring location of these different measures for
% bme performance

if nargin < 1, forceddist = 0; end
if nargin < 1, stationstr = 'stations'; end

strwhat = { 'krig' ; 'RAMP' ; 'diff' ; 'reldiff' };

% getting color bar range
load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','RAMP',floor(forceddist./1000)));
if strcmp(stationstr,'stations'), valplot = stations; else valplot = yearstations; end
for i = 1:length(valplot)
    caxs{i} = [prctile(valplot{i}(:),5) prctile(valplot{i}(:),90)];
end

for l = 1:length(strwhat)
    
    % load file    
    if strcmp(strwhat{l},'diff')
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','RAMP',floor(forceddist./1000)));
        if strcmp(stationstr,'stations'), valplot1 = stations; else valplot1 = yearstations;end  
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','krig',floor(forceddist./1000)));
        if strcmp(stationstr,'stations'), valplot2 = stations; else valplot2 = yearstations;end 
        valplot = cellfun(@(x,y) x-y,valplot1,valplot2,'UniformOutput', false);
    elseif strcmp(strwhat{l},'reldiff')
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','RAMP',floor(forceddist./1000)));
        if strcmp(stationstr,'stations'), valplot1 = stations; else valplot1 = yearstations;end  
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','krig',floor(forceddist./1000)));
        if strcmp(stationstr,'stations'), valplot2 = stations; else valplot2 = yearstations;end 
        valplot = cellfun(@(x,y) (x-y)./y,valplot1,valplot2,'UniformOutput', false);
    else
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat',strwhat{l},floor(forceddist./1000)));
        if strcmp(stationstr,'stations'), valplot = stations; else valplot = yearstations; end        
    end

    % measure names/values
    strnm = {'num of paired modeled and obs';'mean obs';'mean mod';'mean error'; ...
        'mean normalized error';'normalized mean error';'fractional error';'mean absolute error'; ...
        'mean normalized absolute error';'normalized mean absolute error';'fractional absolute error'; ...
        'correlation';'correlation squared';'standard error';'mean squared error'; ...
        'root mean squared error';'normalized root mean squared error'; ...
        'mean error DIV standard error';'mean error squared DIV mean squared error'; ...
        'variance of errors DIV mean squared error';'beta1';'variance of  obs';'varianced of mod'};

    %%% maps overall
    if strcmp(stationstr,'stations') == 1
        for i = 1:length(valplot)

            figure; hold on;

            % country outline
            cd ../09_mfiles_projections
            load('USAcontiguous.mat');
            plotax = ell2lambertcc([x,y],'whiproj2001');
            cd ../14_mfiles_gathervalidation

            % setting axis
            xlabel('km');
            ylabel('km');
            axis([ -3000000 3000000 -2000000 1500000 ]);

            % overlaying the states
            load('../09_mfiles_projections/USAstates5.mat');
            for j = 1:length(X)
                cd ../09_mfiles_projections
                states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
                cd ../14_mfiles_gathervalidation
                plot(states(:,1),states(:,2),'k-');
            end

            % colorplot
            Property={'Marker','MarkerSize','MarkerEdgeColor'};
            Value ={'o',5,[0 0 0]};   
            if i == length(valplot) % deal with beta1
                cax = [prctile(valplot{i},5) prctile(valplot{i},90)];
            elseif length(unique(valplot{i})) == 1            
                cax = [unique(valplot{i})-1 unique(valplot{i})+1];
                valplot{i}(1) = valplot{i}(1) + 0.01;
            elseif l == 1 | l == 2
                cax = caxs{i};
            else
                cax = [prctile(valplot{i},5) prctile(valplot{i},90)];
            end
            colorplot(uniMS,valplot{i},'hot',Property,Value,cax);
            caxis(cax);
            colorbar;
            title(sprintf('%s %d km %s across all years',strwhat{l},floor(forceddist./1000),strnm{i})); 

            % save figure
            set(gcf,'Position',[0 0 800 600]);
            set(gcf,'PaperUnits','inches');    
            set(gcf,'PaperPosition',[0 0 800 600]./100);
            set(gcf,'PaperPositionMode','manual');
            set(gca,'XTickLabel',get(gca,'XTick')/1000);
            set(gca,'YTickLabel',get(gca,'YTick')/1000);
            print(gcf,'-painters','-dpng','-r600',sprintf('figures/maps/%s_%dkm_%s.png', ...
                strwhat{l},floor(forceddist./1000),strnm{i}));

        end
        close all;
    end
    
    %%% maps by year
    if strcmp(stationstr,'yearstations') == 1
        for i = 1:length(valplot)
            for u = 1:length(uniyr)
                disp([l i u]);
                figure; hold on;

                % country outline
                cd ../09_mfiles_projections
                load('USAcontiguous.mat');
                plotax = ell2lambertcc([x,y],'whiproj2001');
                cd ../14_mfiles_gathervalidation

                % setting axis
                xlabel('km');
                ylabel('km');
                axis([ -3000000 3000000 -2000000 1500000 ]);

                % overlaying the states
                load('../09_mfiles_projections/USAstates5.mat');
                for j = 1:length(X)
                    cd ../09_mfiles_projections
                    states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
                    cd ../14_mfiles_gathervalidation
                    plot(states(:,1),states(:,2),'k-');
                end

                % colorplot
                Property={'Marker','MarkerSize','MarkerEdgeColor'};
                Value ={'o',5,[0 0 0]};   
                temp = valplot{i}(u,:);
                idx = ~isnan(temp);
                if length(unique(temp(idx))) == 1            
                    cax = [unique(valplot{i}(u,idx))-1 unique(valplot{i}(u,idx))+1];
                    valplot{i}(u,1) = valplot{i}(u,1) + 0.01;
                elseif l == 1 | l == 2
                    cax = caxs{i};
                else
                    cax = [prctile(valplot{i}(u,:),5) prctile(valplot{i}(u,:),90)];
                    if cax(1)==cax(2), cax(2)=cax(2)+1; end
                end
                colorplot(uniMS(idx,:),valplot{i}(u,idx)','hot',Property,Value,cax);
                caxis(cax);
                colorbar;
                title(sprintf('%s %d km %s for %d',strwhat{l},floor(forceddist./1000),strnm{i},uniyr(u))); 

                % save figure
                set(gcf,'Position',[0 0 800 600]);
                set(gcf,'PaperUnits','inches');    
                set(gcf,'PaperPosition',[0 0 800 600]./100);
                set(gcf,'PaperPositionMode','manual');
                set(gca,'XTickLabel',get(gca,'XTick')/1000);
                set(gca,'YTickLabel',get(gca,'YTick')/1000);
                print(gcf,'-painters','-dpng','-r600',sprintf('figures/maps/%s_%dkm_%d_%s.png', ...
                    strwhat{l},floor(forceddist./1000),uniyr(u),strnm{i}));
            
            end
            close all;
        end
    end

end

end