function [] = analysis()
% this function will call all files in this folder. this folder will 
% gather all of the cross validation results from all the methods and make 
% sense of the following ...
%
% METHODS: kriging, NBDS, StAMP, RAMP
%
% XVAL: LOO, increasing radius, 10 fold, true 10 fold
%
% COMPARISONS: relative difference, absolute difference
%
% DISPLAYS: tables, maps, plots (over time/increasing radius)
%
% METRICS: 1) number of paired modeled and obs, 2) mean obs value,
% 3) mean modeled value, 4) mean bias, 5) normalized bias,
% 6) normalized mean bias, 7) fractional bias, 8) mean error,
% 9) normalized error, 10) normalized mean error, 11) fractional error,
% 12) correlation, 13) correlation squared, 14) standard bias,
% 15) mean squared bias, 16) root mean squared bias,
% 17) normalized root mean squared bias, 18) mean bias/standard bias,
% 19) mean bias squared/mean squared bias, 20) variance of bias/mean squared bias,
% 21) beta1, 22) variance of obs, 23) variance of mod

% loading BME function
cd ../BMELIB2.0b
startup
cd ../14_mfiles_gathervalidation

xstr = { 'LOO' ; '10fold' ; 'true10fold' };
methstr = { 'krig' ; 'RAMP' ; 'CAMP'; 'staticDS' ; 'stDS_add_ind_muli_ind' ; 'stDS_add_dyn_muli_ind' };
radbuff = 0:100000:900000;

% % loo
% for j = [0:1 3:5] % 0:5 % krig, RAMP, CAMP, staticDS, stDS_add_ind_muli_ind, stDS_add_dyn_muli_ind
%     for k = 1:10 % radial buffer % 1:10
%         % method
%         tic
%         [overall,years,regions,yearregions,stations,yearstations, ...
%             nears,yearnears,folds,yearfolds,uniyr,uniR,uniMS,unifold] = stat_MP_bme(j,radbuff(k),0); 
%         save(sprintf('matfiles/allXval_%s_%s_dist%d.mat',xstr{0+1},methstr{j+1},floor(radbuff(k)./1000)), ...
%             'overall','years','regions','yearregions','stations','yearstations', ...
%             'nears','yearnears','folds','yearfolds','uniyr','uniR','uniMS','unifold');
%         toc
%         % creating table of results
%         displayTable(j,radbuff(k),0);
%     end
% end

% 10fold, true 10 fold
for i = 1:2 % 10fold, true 10fold 
    for j = 0:1 % krig, RAMP, 
        % method
        tic
        [overall,years,regions,yearregions,stations,yearstations, ...
            nears,yearnears,folds,yearfolds,uniyr,uniR,uniMS,unifold] = stat_MP_bme(j,radbuff(1),i); 
        save(sprintf('matfiles/allXval_%s_%s_dist%s.mat',xstr{i+1},methstr{j+1},'NA'), ...
            'overall','years','regions','yearregions','stations','yearstations', ...
            'nears','yearnears','folds','yearfolds','uniyr','uniR','uniMS','unifold');
        toc
        % creating table of results
        displayTable(j,radbuff(1),i);
    end
end

% plot figures
% maps
for i = 1:10
    map_MP_bme(radbuff(i),'stations');
    map_MP_bme(radbuff(i),'yearstations');
end
% krig and RAMP statistics by increasing radius/year
map_radius();
for i = 1:10
    map_year(radbuff(i));
end

% add statistics like RMSS?

end