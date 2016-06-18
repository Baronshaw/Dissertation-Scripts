function [RMSE] = myErrFunC(betas,rLags,Crtests,problems)
% this function assumes that A11, A12, and A22 don't change as function of 
% modeled value

n = 10;
lambda = 1/55;

for i = 1:n
    for j = i:n
        
        A11 = betas(1); A12 = betas(2); A22 = betas(3); phi0 = betas(4); phi1 = betas(5);
        xB1 = problems(i); xB2 = problems(j); tau2 = problems(end);
        f{i,j} = (A11-A12.*xB1).*(A11-A12.*xB2).*exp(-1./phi0.*rLags{i,j}) + ...
            (A12-A22.*xB1).*(A12-A22.*xB2).*exp(-1./phi1.*rLags{i,j}) + tau2.*exp(-0.5.*rLags{i,j});
        errors{i,j} = lambda.*sqrt(mean((Crtests{i,j}-f{i,j}').^2));
    
    end
end

temp = errors(:);
temp = cell2mat(temp); 
RMSE = sum(temp);

end