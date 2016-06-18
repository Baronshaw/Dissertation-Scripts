function [] = displayXvalStatistics_10fold(whichXval)
% this function will display the results of the 10 fold cross validation
% statistics

if nargin < 1, whichXval = 1; end

% 1 = ME, 2 = MAE, 3 = MSE, 4 = pear, 5 = spear

switch whichXval
    case 1, strXval = 'ME';
    case 2, strXval = 'MAE';
    case 3, strXval = 'MSE';
    case 4, strXval = 'r^{2}';
    case 5, strXval = 'spear r^{2}';
end

% displaying results in chart form
methods = [ 0 1 1 ; 0 0 1 ; 1 0 1 ];
[r c] = size(methods);
n = 1;
for i = 1:10

    folds(n,1) = i;

    for j = 1:r            
        % loading results
        if methods(j,1) == 0, softstr = '_nosoft'; else softstr = '_soft'; end
        if methods(j,2) == 1, constr = '_constant'; else constr = '_long'; end
        if methods(j,3) == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end     
        load(sprintf('../matfiles/Xval_10fold%s%s%s_results.mat',softstr,constr,gaussstr));
        if whichXval == 4 | whichXval == 5 
            Xvalstat(n,j) = (alltogether{whichXval}).^2; 
        else
            Xvalstat(n,j) =alltogether{whichXval}; 
        end
    end 

    n = n + 1;  
        
end

% outid = fopen(sprintf('Xval_results/showXval_%s.csv',strXval),'w+');
% fprintf(outid,'%s');
% fclose(outid);
dlmwrite(sprintf('showXval_%s.csv',strXval),[folds Xvalstat],'delimiter',',','precision',6,'-append','roffset',1)    

end