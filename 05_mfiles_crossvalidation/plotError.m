function [] = plotError(soft,constant,gauss,daydisp)
% this function will plot figures of error for a given day across
% the US

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not
if nargin < 4, daydisp = [2001 7 1]; end

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% load cross-valdiation results
load(sprintf('Xval_10fold%s%s%s_results_test1.mat',softstr,constr,gaussstr));

% variables to work with
unisID = unique(ckall(:,1:2),'rows');
errall = zkall-zhall;
ckdisp = datevec(ckall(:,3));
idx = daydisp(1) == ckdisp(:,1) & daydisp(2) == ckdisp(:,2) & daydisp(3) == ckdisp(:,3);
dayshow = daydisp(1)*10000 + daydisp(2)*100 + daydisp(3);

figure; hold on;

% country outline
cd ../09_mfiles_projections
load('USAcontiguous.mat');
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../05_mfiles_crossvalidation

% setting axis
xlabel('km');
ylabel('km');
axis([ -3000000 3000000 -2000000 1500000 ]);

% overlaying the states
load('../09_mfiles_projections/USAstates5.mat');
for j = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    plot(states(:,1),states(:,2),'k-');
end

% colorplot
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'o',5,[0 0 0]};
colorplot(ckall(idx,1:2),errall(idx),'hot',Property,Value);
cax = [prctile(errall,5) prctile(errall,95)];
caxis(cax);
colorbar;
title(sprintf('Error on %d %s %s %s %s',dayshow,softstr(2:end),constr(2:end),gaussstr(2:end))); 

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
print(gcf,'-painters','-dpdf','-r600',sprintf('abserror%d_%s%s%s_test1.pdf',dayshow,softstr,constr,gaussstr));
print(gcf,'-painters','-dpng','-r600',sprintf('abserror%d_%s%s%s_test1.png',dayshow,softstr,constr,gaussstr));

end