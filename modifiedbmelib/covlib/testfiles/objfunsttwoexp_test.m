function objfunsttwoexp_test()

sptllag = 0:15;
sptlcovval = exponentialC(sptllag,[0.7,4]) + ...
             exponentialC(sptllag,[0.3,8]);
sptlweight  = ones(size(sptllag)) + randn(size(sptllag)) * 0.1;

templag = 0:10:100;
tempcovval = exponentialC(templag,[0.7,40]) + ...
             exponentialC(templag,[0.3,80]);
tempweight  = ones(size(templag)) + randn(size(templag)) * 0.1;

covparam1 = 0.1:0.1:0.9; % 9 elements
covparam2 = 2:10;        % 9 elements
covparam3 = 20:10:100;   % 9 elements
covparam4 = 0.2:0.1:0.4; % 3 elements
covparam5 = 2:10;        % 9 elements
covparam6 = 20:10:100;   % 9 elements

fval = ones(size(covparam1,2) * size(covparam2,2) *...
            size(covparam3,2) * size(covparam4,2) *...
            size(covparam5,2)*size(covparam6,2),1) * NaN;
covvals = [];

midx = 1;
for i = 1:size(covparam1,2)
    disp(['i = ',num2str(i),'/9']);
    for j = 1:size(covparam2,2)
        for k = 1:size(covparam3,2)
            for l = 1:size(covparam4,2)
                for m = 1:size(covparam5,2)
                    for n = 1:size(covparam6,2)
                        covparam = [covparam1(i),covparam2(j),...
                                    covparam3(k),covparam4(l),...
                                    covparam5(m),covparam6(n)];
                        fval(midx) = objfunsttwoexp(covparam,sptllag,...
                                         sptlcovval,sptlweight,templag,...
                                         tempcovval,tempweight);
                        covvals = [covvals;covparam];
                        midx = midx + 1;
                    end
                end
            end
        end
    end
end

disp(['index of min(f): ',num2str(transpose(find(fval==min(fval))))]);
disp(['min(f)         : ',num2str(min(fval))]);

disp(['covparam1 (Partial Sill1): ',...
      num2str(transpose(covvals(fval == min(fval),1)))]);
disp(['covparam2 (Sptl Range1)  : ',...
      num2str(transpose(covvals(fval == min(fval),2)))]);
disp(['covparam3 (Temp Range1)  : ',...
      num2str(transpose(covvals(fval == min(fval),3)))]);
disp(['covparam4 (Partial Sill2): ',...
      num2str(transpose(covvals(fval == min(fval),4)))]);
disp(['covparam5 (Sptl Range2)  : ',...
      num2str(transpose(covvals(fval == min(fval),5)))]);
disp(['covparam6 (Temp Range2)  : ',...
      num2str(transpose(covvals(fval == min(fval),6)))]);
