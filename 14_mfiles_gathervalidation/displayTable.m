function [] = displayTable(soft,forceddist,loo)
% this function will create data tables for the results

% parameters
if nargin < 1, soft = 1; end 
if nargin < 2, forceddist = 0; end
if nargin < 3, loo = 1; end

if soft == 0, softstr = 'krig'; 
elseif soft == 1, softstr = 'RAMP'; 
elseif soft == 2, softstr = 'CAMP'; 
elseif soft == 3, softstr = 'staticDS'; 
elseif soft == 4, softstr = 'stDS_add_ind_muli_ind'; 
elseif soft == 5, softstr = 'stDS_add_dyn_muli_ind';   
end
if loo == 0, xstr = 'LOO'; elseif loo == 1, xstr = '10fold'; else xstr = 'true10fold'; end

% loading matfile
if loo == 0
    load(sprintf('matfiles/allXval_%s_%s_dist%d.mat',xstr,softstr,floor(forceddist./1000)));
elseif loo == 1 | loo == 2
    load(sprintf('matfiles/allXval_%s_%s_dist%s.mat',xstr,softstr,'NA'));
end

% header string
str1 = 'num pairs,mean obs,mean mod,mean bias,norm bias,norm mean bias,';
str2 = 'frac bias,mean err,norm err,norm mean err,frac err,corr,corr squared,';
str3 = 's bias,ms bias,rms bias,nrms bias,mean bias DIV s bias,mean bias squared DIV ms bias,v bias DIV ms bias,beta1,v obs,v mod';
strval = sprintf('%s%s%s',str1,str2,str3);

% overall
outid = fopen(sprintf('tables/overall_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('tables/overall_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)), ...
    overall,'delimiter',',','precision',6,'-append','roffset',1)

% years
outid = fopen(sprintf('tables/years_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)),'w+');
fprintf(outid,'years,%s',strval);
fclose(outid);
dlmwrite(sprintf('tables/years_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)), ...
    [uniyr years],'delimiter',',','precision',6,'-append','roffset',1)

% regions
outid = fopen(sprintf('tables/regions_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)),'w+');
fprintf(outid,'regions,%s',strval);
fclose(outid);
dlmwrite(sprintf('tables/regions_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)), ...
    [uniR regions],'delimiter',',','precision',6,'-append','roffset',1)

% yearregions
allyr = repmat(uniyr',length(uniR),1);
allR = repmat(uniR,length(uniyr),1);
allYR = cellfun(@(x) x(:), yearregions, 'UniformOutput', false);
outid = fopen(sprintf('tables/yearregions_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)),'w+');
fprintf(outid,'years,regions,%s',strval);
fclose(outid);
dlmwrite(sprintf('tables/yearregions_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)), ...
    [allyr(:) allR allYR],'delimiter',',','precision',6,'-append','roffset',1)
    
% nears
allN = [1:10]';
outid = fopen(sprintf('tables/nears_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)),'w+');
fprintf(outid,'nears,%s',strval);
fclose(outid);
dlmwrite(sprintf('tables/nears_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)), ...
    [allN nears],'delimiter',',','precision',6,'-append','roffset',1)

% yearnears
allyr = repmat(uniyr',10,1);
allN = repmat([1:10]',length(uniyr),1);
allYN = cellfun(@(x) x(:), yearnears, 'UniformOutput', false);
outid = fopen(sprintf('tables/yearnears_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)),'w+');
fprintf(outid,'years,nears,%s',strval);
fclose(outid);
dlmwrite(sprintf('tables/yearnears_%s_%s_dist%d.csv',xstr,softstr,floor(forceddist./1000)), ...
    [allyr(:) allN allYN],'delimiter',',','precision',6,'-append','roffset',1)

% folds
if loo == 1 | loo == 2
    outid = fopen(sprintf('tables/folds_%s_%s_dist%d.csv',xstr,softstr,0),'w+');
    fprintf(outid,'fold,%s',strval);
    fclose(outid);
    dlmwrite(sprintf('tables/folds_%s_%s_dist%d.csv',xstr,softstr,0), ...
        [unifold folds],'delimiter',',','precision',6,'-append','roffset',1)
end

% yearfolds
if loo == 1 | loo == 2
    allyr = repmat(uniyr',10,1);
    allF = repmat([1:10]',length(uniyr),1);
    allYF = cellfun(@(x) x(:), yearfolds, 'UniformOutput', false);
    outid = fopen(sprintf('tables/yearfolds_%s_%s_dist%d.csv',xstr,softstr,0),'w+');
    fprintf(outid,'years,folds,%s',strval);
    fclose(outid);
    dlmwrite(sprintf('tables/yearfolds_%s_%s_dist%d.csv',xstr,softstr,0), ...
        [allyr(:) allF allYF],'delimiter',',','precision',6,'-append','roffset',1)
end

end