function [] = getMeanTrend(years)
% this function will calculate the mean trend and save the results. Note:
% only 2001 data was used to calculate the mean trend

if nargin < 1, years = 2001; end

% loading 2001 data
load(sprintf('../matfiles/prepObs_%d.mat',years));
yrs = floor(yrmodaObs./10000);
mos = floor( (yrmodaObs-yrs*10000) ./ 100 );
das = yrmodaObs - yrs*10000 - mos*100;
pd = [coordObs datenum(yrs,mos,das)];
zd = Obs;

% find the maximum of closest neighbors
[Zd,MSd,MEd,nanratio]=valstv2stg(pd,zd); MSd = [MSd(:,2) MSd(:,1)];
DMS = sqrt(bsxfun(@plus,dot(MSd,MSd,2),dot(MSd,MSd,2)')-2*(MSd*MSd'));
DME = abs(bsxfun(@minus,MEd,MEd'))';
temp = sort(DMS,2);
maxMinNei = max(temp(:,2));

% all the smoothing parameters to try
ars = [ 20 40 20000 40000 50000:50000:2500000 ];
ats = [ 4:4:20 25:25:300 ];

spatRad = max([2*ars' maxMinNei*ones(length(ars),1)],[],2);
atconstant = 30;
arconstant = 300000;
smoothingParamsSpace = [spatRad ars' 2*atconstant*ones(length(ars),1) atconstant*ones(length(ars),1)];
smoothingParamsTime = [spatRad(1)*ones(length(ats),1) arconstant*ones(length(ats),1) 2*ats' ats'];

% because this is the testing portion of the mean trend, I don't need to
% calculate the mean trend for ALL the data, just for some days and some
% locations so I can visualize the results in order to determine the
% optimal parameters

% for space, I'll look at one day each month. First see which day in each
% month measured the most
temp = histc(pd(:,3),[datenum(years,1,1):datenum(years,12,31)]');
daysinyear = datevec([datenum(years,1,1):datenum(years,12,31)]');
for i = 1:12
    idx = daysinyear(:,2) == i;
    dummy = find(temp(idx)==max(temp(idx)));
    daysinmonth = daysinyear(idx,:);
    maxshow(i,:) = daysinmonth(dummy(1),1:3);
end
[lia lib] = ismember(pd(:,3),datenum(maxshow)); 
pI = pd(lia,:);

for i = 1:size(smoothingParamsSpace,1)
    smoothingParam = smoothingParamsSpace(i,:);
    cd ../10_mfiles_newmeantrend
    tic
    [mI]=expKernelSmooth_stv(pd,zd,smoothingParam,pI);
    toc
    cd ../02_mfiles_offset  
    % saving results
    save(sprintf('../matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',years,smoothingParam), ...
        'pd','zd','smoothingParam','pI','mI');
end

% for time, I'll look at 10 locations with the most amount of data
uniloc = unique(pd(:,1:2),'rows');
moncnt = NaN*ones(size(uniloc,1),1);
for i = 1:size(uniloc,1)
    idx = uniloc(i,1) == pd(:,1) & uniloc(i,2) == pd(:,2);
    moncnt(i) = sum(idx);
end
[sortm sortidx] = sort(moncnt,'descend');
[lia lib] = ismember(pd(:,1:2),uniloc(sortidx(1:10),1:2),'rows'); 
pI = pd(lia,:);

for i = 1:size(smoothingParamsTime,1)
    smoothingParam = smoothingParamsTime(i,:);
    cd ../10_mfiles_newmeantrend
    tic
    [mI]=expKernelSmooth_stv(pd,zd,smoothingParam,pI);
    toc
    cd ../02_mfiles_offset   
    % saving results
    save(sprintf('../matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',years,smoothingParam), ...
        'pd','zd','smoothingParam','pI','mI');
end

% % final smoothing parameters
% ar2 = [20000 50000 300000 1000000];
% at2 = [10 20 50 200];

end