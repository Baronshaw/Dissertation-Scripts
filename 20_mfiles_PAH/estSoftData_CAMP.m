function [] = estSoftData_CAMP(pah)
% estimate soft data using an pseudo one-the-fly CAMP method

if nargin < 1, pah = 1; end

% load data
load('matfiles/pah_data.mat');

% defining variables
idx = ~isnan(val(:,pah+4));
zPAH = log(val(idx,pah+4));
zPMsub = log(val(idx,4));
zPM = log(val(:,4));
cPM = [ProjectX ProjectY Time];
cPAH = [ProjectX(idx) ProjectY(idx) Time(idx)];
buff = std(zPMsub);

% mean and variance for the soft data locations
mPAHs = NaN*ones(length(cPM),1);
mPAHs(idx) = zPAH;
vPAHs = NaN*ones(length(cPM),1);
vPAHs(idx) = 0;
for i = 1:length(cPM)
    
    if isnan(mPAHs(i)) % if soft data needs to be calculated
        
        idx = zPMsub >= zPM(i)-buff & zPMsub <= zPM(i)+buff;

        n = 2;
        while sum(idx) < 2
            buffmod = n.*buff;
            idx = zPMsub >= zPM(i)-buffmod & zPMsub <= zPM(i)+buffmod;
            n = n + 1;
        end

        mPAHs(i) = mean(zPAH(idx));
        vPAHs(i) = var(zPAH(idx));
        
    end
    
end

% saving soft data
save(sprintf('matfiles/soft_%s_CAMP.mat',valname{pah+4}), ...
    'cPM','mPAHs','vPAHs','cPAH','zPAH','zPM','zPMsub','buff'); 

end