function [RMSE] = myErrFunCST(betas,rLags,Crtests,problems)
% this function assumes that A11, A12, and A22 don't change as function of 
% modeled value for the covarinace model of the space/time downscaler

n = 10;
lambda = 1/55;
T = length(rLags);
f = cell(T,1);
errors = cell(T,1);

for i = 1:T
    
    f{i} = cell(n,n);
    errors{i} = cell(n,n);
    
    for j = 1:n
        for k = j:n

            A11 = betas(1); A12 = betas(2); A22 = betas(3); phi0 = betas(4); phi1 = betas(5);
            xB1 = problems{i}(j); xB2 = problems{i}(k); tau2 = problems{i}(end);
            f{i}{j,k} = (A11-A12.*xB1).*(A11-A12.*xB2).*exp(-1./phi0.*rLags{i}{j,k}) + ...
                (A12-A22.*xB1).*(A12-A22.*xB2).*exp(-1./phi1.*rLags{i}{j,k}) + tau2.*exp(-0.5.*rLags{i}{j,k});
            if ~isempty(Crtests{i}{j,k})
                errors{i}{j,k} = lambda.*sqrt(mean((Crtests{i}{j,k}-f{i}{j,k}').^2));
            end

        end
    end
end

temp = cellfun(@(x) cell2mat(x(:)),errors,'UniformOutput',false); 
temp = cell2mat(temp(:));
RMSE = sum(temp);

end