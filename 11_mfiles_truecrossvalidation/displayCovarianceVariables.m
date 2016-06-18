function [] = displayCovarianceVariables()
% this function will display the values of the covariance parameters in a
% .csv file

matall = NaN*ones(11,6); % these are hard-coded
matall(:,1) = 0:10; % naming the folds

% loading old parameters
load('../matfiles/covmod_r_long_joint exponential exponential_joint.mat');
matall(1,2:end) = [f.alp f.ar1 f.ar2 f.at1 f.at2];   

% getting each fold
for i = 1:10
    
    load(sprintf('../matfiles/covmod_true10fold_%d.mat',i));
    matall(i+1,2:end) = [f.alp f.ar1 f.ar2 f.at1 f.at2];
   
end

matall(:,3:4) = floor(matall(:,3:4)./1000);

% writing file
strval = 'fold,alpha,ar1,ar2,at1,at2';
outid = fopen(sprintf('covariance_values_true10fold.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('covariance_values_true10fold.csv'), ...
    matall,'delimiter',',','precision',6,'-append','roffset',1)

end