function [] = spacetimeDSXval(additive,multiplicative)
% this function will perform the LOOCE method for the space/time downscaler

% only use the following if parallel computing is needed
% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "spacetimeDSXval" -logfile "runall_Cluster12.out"

if nargin < 1, additive = 'ind'; end 
if nargin < 2, multiplicative = 'ind'; end

% loading BME function
cd ../BMELIB2.0b
startup
cd ../12_mfiles_othermethods

% load parameter files
load(sprintf('matdata/STDS_results_add_%s_muli_%s.mat',additive,multiplicative));

% kriging parameters
cs = []; zs = []; vs = [];
nsmax = 0;
order = NaN;
options = BMEoptions;
options(3) = 150000;

% looping through each station
len = length(unique(ch,'rows'));
unich = unique(ch,'rows');
zk_mtr = cell(len,1); % estimated mean
vk = cell(len,1); % estimated variance
zk = cell(len,1); % estimated mean with mean trend added back in
zhX = cell(len,1); % observed values
zhX_mtr = cell(len,1); % observed values with mean trend removed
ck = cell(len,1); % locations of X val
for i = 1:len
    disp([i len]);
 
    idx = unich(i,1) == ch(:,1) & unich(i,2) == ch(:,2);
    zhX{i} = zho(idx);
    zhX_mtr{i} = zho_mtr(idx);
    ck{i} = [ch(idx,:) cht(idx);];

    covmodel = {'nuggetC/nuggetC','exponentialC_DS/nuggetC','exponentialC_DS/nuggetC'};
    zhmp = prctile(zhm,0:10:100);
    problems = NaN*ones(1,10);
    for k = 1:10
        problems(k) = (zhmp(k+1)-zhmp(k))./2 + zhmp(k); 
    end
    xB1 = problems(1); xB2 = problems(2); 
    covparam = {[tau2 tau2] [(A11-A12.*xB1).*(A11-A12.*xB2) 3.*phi0 tau2] [(A12-A22.*xB1).*(A12-A22.*xB2) 3*phi1 tau2]};
    dmax3 = (3.*phi0+3.*phi1)/(tau2+tau2);
    dmax = [2000000 20 dmax3];

    % find modeled values of the observed data
    temp1 = arrayfun(@(x) find(abs(x-zhmp)==min(abs(x-zhmp))),zhm(idx),'UniformOutput',false);
    temp2 = find( cell2mat(cellfun(@(x) length(x)>1,temp1,'UniformOutput',false)) == 1);
    for k = 1:length(temp2)
        temp1{temp2(k)} = temp1{temp2(k)}(1);
    end
    zhmC = zhmp(cell2mat(temp1));
    xB1M = repmat(zhmC,length(zhmC),1); xB2M = repmat(zhmC',1,length(zhmC));
    covparamd_DS = {[tau2] [(A11-A12.*xB1M).*(A11-A12.*xB2M)] [(A12-A22.*xB1M).*(A12-A22.*xB2M)]};

    % find modeled values of the estimation locations
    if ismember(ck{i},[ch cht],'rows') == 1 % estimation location in the data 
        [ia ib] = ismember(ck{i},[ch cht],'rows');
    else % new estimation location
        disp('Can''t account for new location in crossvalidation');
    end
    temp1 = arrayfun(@(x) find(abs(x-zhmp)==min(abs(x-zhmp))),zhm(ib),'UniformOutput',false);
    temp2 = find( cell2mat(cellfun(@(x) length(x)>1,temp1,'UniformOutput',false)) == 1);
    for k = 1:length(temp2)
        temp1{temp2(k)} = temp1{temp2(k)}(1);
    end
    zkm = zhmp(cell2mat(temp1));
    covparamk_DS = {[tau2] [(A11-A12.*zkm).*(A11-A12.*zhmC)] [(A12-A22.*zkm).*(A12-A22.*zhmC)]};
    covparamk0_DS = {[tau2] [(A11-A12.*zkm).*(A11-A12.*zkm)] [(A12-A22.*zkm).*(A12-A22.*zkm)]};

    % calculate an estimate
    nhmax = 70; 
    if sum(idx) < nhmax, nhmax = sum(idx); end
    chX = [ch(~idx,:) cht(~idx)];
    zh_mtr = zho_mtr(~idx);
    [zk_mtr{i},vk{i}]=krigingME_DS(ck{i},chX,cs,zh_mtr,zs,vs,covmodel,covparam, ...
        covparamd_DS,covparamk_DS,covparamk0_DS,nhmax,nsmax,dmax,order,options);
    zk{i} = zk_mtr{i} + zkm'; % add mean trend 

end

% saving LOOCE
save(sprintf('matdata/LOOCE_DS_add_%s_muli_%s.mat',additive,multiplicative), ...
    'unich','zk_mtr','vk','zk', 'zhX', 'zhX_mtr','ck');

end