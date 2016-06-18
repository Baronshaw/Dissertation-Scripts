function [] = tableSimu()
% this function will create a table showing the correlation coefficient
% across all grids for July 1, 2001 between lambda1 and lambda1*, 
% lambda2 and lambda2*, ME and ME*, and SE and SE* 
% across the Constant, CAMP, and RAMP method 

% selected lambda1 and lambda2
load('PM2p5_meanGivMod_yr2001.mat');

% pick day
unidays = datenum(2001,7,1);
idx = css(:,3) == unidays;

% CONSTANT - recalculated lambda1 and lambda2
load('recalculatelambdaGivMod_Constant.mat');
% lambda1_recalcModGrid, lambda2_recalcsimuModGrid

col1 = [ corr(meanGivCMAQ_all(idx),real(lambda1_recalcsimuModGrid(idx))) ; 
    corr(varGivCMAQ_all(idx),lambda2_recalcsimuModGrid(idx)) ;
    corr(dailyCTMvorder(idx)-meanGivCMAQ_all(idx),dailyCTMvorder(idx)-real(lambda1_recalcsimuModGrid(idx))) ];

% CAMP - recalculated lambda1 and lambda2
load('recalculatelambdaGivMod_CAMP.mat');
% lambda1_recalcModGrid, lambda2_recalcsimuModGrid

col2 = [ corr(meanGivCMAQ_all(idx),real(lambda1_recalcsimuModGrid(idx))) ; 
    corr(varGivCMAQ_all(idx),lambda2_recalcsimuModGrid(idx)) ;
    corr(dailyCTMvorder(idx)-meanGivCMAQ_all(idx),dailyCTMvorder(idx)-real(lambda1_recalcsimuModGrid(idx))) ];

% CAMP - recalculated lambda1 and lambda2
load('recalculatelambdaGivMod.mat');
% lambda1_recalcModGrid, lambda2_recalcsimuModGrid

col3 = [ corr(meanGivCMAQ_all(idx),real(lambda1_recalcsimuModGrid(idx))) ; 
    corr(varGivCMAQ_all(idx),lambda2_recalcsimuModGrid(idx)) ;
    corr(dailyCTMvorder(idx)-meanGivCMAQ_all(idx),dailyCTMvorder(idx)-real(lambda1_recalcsimuModGrid(idx))) ];

% column string
colstr = {'lambda1','lambda2','ME=CMAQ-lambda1'};

% create table
fileID = fopen('tables/correlation_simulation.csv','w+');
fprintf(fileID,'%s,%s,%s,%s\n','parameter','Constant','CAMP','RAMP');
for i = 1:length(colstr)
    fprintf(fileID,'%s,%f,%f,%f\n',colstr{i},col1(i),col2(i),col3(i));
end
fclose(fileID);

end