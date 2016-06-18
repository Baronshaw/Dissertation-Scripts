function [] = displayXvalforcedisolation(years,soft,constant,gauss,percent,forceddist)

if nargin < 1, years = 2001; end % year of time series
if nargin < 2, soft = 1; end % soft data or not
if nargin < 3, constant = 0; end % constant offset or not
if nargin < 4, gauss = 1; end % gaussian soft data or not
if nargin < 5, percent = 2; end
if nargin < 6, foreddist = 200; end % distance in kilometers

if soft == 0, softstr = '_nosoft'; else softstr = ''; end
if constant == 1, constr = '_constant'; else constr = ''; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = ''; end

load(sprintf('Xval_results/XvalforcedisolationTESTING_per%dfor%d_%d%s%s%s.mat', ...
    percent,foreddist,years,softstr,constr,gaussstr));
    
SE1 = (zkresoff-subZresoff).^2;

load XvalmajorclustersTESTING_2001.mat
SE = (zkresoff-subZresoff).^2;

idxmin = SE<SE1;

blah = ck(idxmin,:);
uniloc = unique(ck(:,1:2),'rows');


for i = 1:length(uniloc)
    disp1 = sum(uniloc(i,1)==ck(:,1)&uniloc(i,2)==ck(:,2));
    disp2 = sum(uniloc(i,1)==blah(:,1)&uniloc(i,2)==blah(:,2));
    disp3(i) = 100*disp2/disp1;
    % disp(disp3(i));
end

idx = disp3>75; % stations where soft data helped at least 75% of the time

locsoft = uniloc(idx,:);

% maybe plot these stations (?)
% plot all stations
% plot Xval stations
% plot stations where soft data helped in estimation

% picking the cross-validation locations
cd ..
load('distprepsoft_new.mat');

% idx for Xval removal
idxremove = cMSObs{1}(:,1) > -1600000 & cMSObs{1}(:,1) < -900000 & ...
    cMSObs{1}(:,2) > 200000 & cMSObs{1}(:,2) < 650000;

figure; hold on;
plot(cMSObs{1}(:,1),cMSObs{1}(:,2),'bo', ...
    cMSObs{1}(idxremove,1),cMSObs{1}(idxremove,2),'co', ...
    locsoft(:,1),locsoft(:,2),'ro');

% blue = monitoring stations
% red + cyan = removed monitoring stations
% cyan = stations where squared error of scenario 3 is less than squared
% error of scenario 2 at least 75% of the time


end