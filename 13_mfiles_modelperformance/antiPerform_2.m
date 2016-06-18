function [] = antiPerform_2(forceddist)
% this function will test the claims about anticipating performance and
% write results into a table

% What's most important about this script is to find the maximum difference
% between kriging and RAMP when looking at MSE and R2

% ideas: go through all the MP statistics
% combine metrics with FRM/TEOM/East/West/STN/IMPROVE/URBAN/RURAL and MP
% metrics
% different dividing lines (exclude parts of the country)

if nargin < 1, forceddist = 0; end

% loading info
load(sprintf('matfiles/allInfo_%dkm.mat',floor(forceddist./1000)));

% get region
if ~exist('matfiles/inregion.mat')
    allregions = shaperead('FWS_LCC/FWS_LCC.shp');
    for k = 1:length(allregions)
        tic
        cd ../09_mfiles_projections
        allregions_p{k,1} = ell2lambertcc([allregions(k).X',allregions(k).Y'],'whiproj2001');
        cd ../13_mfiles_modelperformance
        in{k,1} = inpolygon(ck_base(:,1),ck_base(:,2),allregions_p{k,1}(:,1),allregions_p{k,1}(:,2));
        toc
    end
    save('matfiles/inregion.mat','in','allregions');
else
    load('matfiles/inregion.mat')
end

% get season
[yr mo da] = datevec(ck_base(:,3));
IsWinter = NaN*ones(length(ck_base),1); IsWinter(mo==1|mo==2|mo==12) = 1;
IsSpring = NaN*ones(length(ck_base),1); IsSpring(mo==3|mo==4|mo==5) = 1;
IsSummer = NaN*ones(length(ck_base),1); IsSummer(mo==6|mo==7|mo==8) = 1;
IsFall = NaN*ones(length(ck_base),1); IsFall(mo==9|mo==10|mo==11) = 1;

% get longitude
cd ../09_mfiles_projections
testing = [-125:5:-90]'; 
M = ell2lambertcc([testing 40*ones(length(testing),1)],'whiproj2001');
cd ../13_mfiles_modelperformance
Is125 = NaN*ones(length(ck_base),1); Is125(ck_base(:,1)<M(1,1)) = 1;
Is120 = NaN*ones(length(ck_base),1); Is120(ck_base(:,1)<M(2,1)) = 1;
Is115 = NaN*ones(length(ck_base),1); Is115(ck_base(:,1)<M(3,1)) = 1;
Is110 = NaN*ones(length(ck_base),1); Is110(ck_base(:,1)<M(4,1)) = 1;
Is105 = NaN*ones(length(ck_base),1); Is105(ck_base(:,1)<M(5,1)) = 1;
Is100 = NaN*ones(length(ck_base),1); Is100(ck_base(:,1)<M(6,1)) = 1;
Is95 = NaN*ones(length(ck_base),1); Is95(ck_base(:,1)<M(7,1)) = 1;
Is90 = NaN*ones(length(ck_base),1); Is90(ck_base(:,1)<M(8,1)) = 1;

% going through every combination of locational info
idx_all = { IsFRM ; IsTEOM ; IsitWest ; ~IsitWest ; IsSTN ; IsIMPROVE ; ...
    IsUrban ; IsRural ; IsSuburban ; ...
    IsWinter ; IsSpring ; IsSummer ; IsFall ; ...
    Is125 ; Is120 ; Is115 ; Is110 ; Is105 ; Is100 ; Is95 ; Is90 ; Dist2NMon > 60000 ; ...
    m2DmsBias_2 >= 0.75 ; m2DmsBias_2 < 0.25 ; ...
    msBias_2 < 10 ; msBias_2 > 100 ; ...
    (mBias_2).^2 < prctile((mBias_2).^2,10) ; (mBias_2).^2 > prctile((mBias_2).^2,90) ; ...
    (sBias_2).^2 < prctile((sBias_2).^2,10) ; (sBias_2).^2 > prctile((sBias_2).^2,90) ; ...
    s2DmsBias_2 >= 0.75 ; s2DmsBias_2 < 0.25 ; ...
    nmErr_2 < prctile(nmErr_2,10) ; nmErr_2 > prctile(nmErr_2,90) ; ...
    vMod_2 < prctile(vMod_2,10) ; vMod_2 > prctile(vMod_2,90) ; 
    in{1}==1 ; in{2}==1; in{3}==1; in{4}==1 ; in{5}==1 ; in{6}==1 ; ...
    in{7}==1 ; in{8}==1 ; in{9}==1 ; in{10}==1 ; in{11}==1 ; in{12}==1 ; ...
    in{13}==1 ; in{14}==1 ; in{15}==1 ; in{16}==1 ; in{17}==1 ; in{18}==1 ; ...
    in{19}==1 ; in{20}==1 ; in{21}==1 ; in{22}==1 ; in{23}==1 ; in{24}==1 };
idx_all_str = { 'IsFRM' ; 'IsTEOM' ; 'IsitWest' ; '~IsitWest' ; 'IsSTN' ; 'IsIMPROVE' ; ...
    'IsUrban' ; 'IsRural' ; 'IsSuburban' ; ...
    'IsWinter' ; 'IsSpring' ; 'IsSummer' ; 'IsFall' ; ...
    'Is125' ; 'Is120' ; 'Is115' ; 'Is110' ; 'Is105' ; 'Is100' ; 'Is95' ; 'Is90' ; 'Dist2NMon > 60000' ; ...
    'm2DmsBias_2 >= 0.75' ; 'm2DmsBias_2 < 0.25' ; ...
    'msBias_2 < 10' ; 'msBias_2 > 100' ; ...
    '(mBias_2).^2 < prctile((mBias_2).^2 10)' ; '(mBias_2).^2 > prctile((mBias_2).^2 90)' ; ...
    '(sBias_2).^2 < prctile((sBias_2).^2 10)' ; '(sBias_2).^2 > prctile((sBias_2).^2 90)' ; ...
    's2DmsBias_2 >= 0.75' ; 's2DmsBias_2 < 0.25' ; ...
    'nmErr_2 < prctile(nmErr_2 10)' ; 'nmErr_2 > prctile(nmErr_2 90)' ; ...
    'vMod_2 < prctile(vMod_2 10)' ; 'vMod_2 > prctile(vMod_2 90)' ; 
    'in{1}==1' ; 'in{2}==1'; 'in{3}==1'; 'in{4}==1' ; 'in{5}==1' ; 'in{6}==1' ; ...
    'in{7}==1' ; 'in{8}==1' ; 'in{9}==1' ; 'in{10}==1' ; 'in{11}==1' ; 'in{12}==1' ; ...
    'in{13}==1' ; 'in{14}==1' ; 'in{15}==1' ; 'in{16}==1' ; 'in{17}==1' ; 'in{18}==1' ; ...
    'in{19}==1' ; 'in{20}==1' ; 'in{21}==1' ; 'in{22}==1' ; 'in{23}==1' ; 'in{24}==1' };
for k = 1:2        
    for i = 1:length(idx_all)  
        idx = idx_all{i} == 1 & ~isnan(Mod_2);
        if k == 1, Mod = zk_hard(idx); else Mod = zk_soft(idx); end
        Obs = zh_base(idx);
        if sum(idx)>0, all_R2(i,k) = (corr(Mod,Obs)).^2; else all_R2(i,k) = NaN; end
        all_msBias(i,k) = mean((Mod-Obs).^2);
        all_counter{i,k} = [i k];
        all_sumidx(i,k) = sum(idx);
    end
end

% overall
fid = fopen(sprintf('tables/antiPerfrom_%dkm.csv',floor(forceddist./1000)),'wt');
fprintf(fid,'names,diff R2,diff MSE,count\n');
for i = 1:length(idx_all_str)
    fprintf(fid,'%s,%f,%f,%f\n',idx_all_str{i}, ...
        all_R2(i,2)-all_R2(i,1),all_msBias(i,2)-all_msBias(i,1),all_sumidx(i,1));
end
fclose(fid);

for k = 1:2
    n = 1;
    for i = 1:length(idx_all)-1
        for j = i+1:length(idx_all)
            idx = idx_all{i} == 1 & idx_all{j} == 1 & ~isnan(Mod_2);
            if k == 1, Mod = zk_hard(idx); else Mod = zk_soft(idx); end
            Obs = zh_base(idx);
            if sum(idx)>0, all_pairs_R2(n,k) = (corr(Mod,Obs)).^2; else all_pairs_R2(n,k) = NaN; end
            all_pairs_msBias(n,k) = mean((Mod-Obs).^2);
            all_pairs_counter{n,k} = [i j k];
            all_pairs_sumidx(n,k) = sum(idx);
            str_col1{n,1} = idx_all_str{i};
            str_col2{n,1} = idx_all_str{j};
            n = n + 1;
        end
    end
end

% pairs
fid = fopen(sprintf('tables/antiPerfrom_pairs_%dkm.csv',floor(forceddist./1000)),'wt');
fprintf(fid,'names1,names2,diff R2,diff MSE,change R2,change MSE,count\n');
for i = 1:n-1
    fprintf(fid,'%s,%s,%f,%f,%f,%f,%f\n',str_col1{i},str_col2{i}, ...
        all_pairs_R2(i,2)-all_pairs_R2(i,1),all_pairs_msBias(i,2)-all_pairs_msBias(i,1), ...
        (all_pairs_R2(i,2)-all_pairs_R2(i,1))./all_pairs_R2(i,1), ...
        (all_pairs_msBias(i,2)-all_pairs_msBias(i,1))./all_pairs_msBias(i,1), ...
        all_pairs_sumidx(i,1));
end
fclose(fid);

testing1 = (all_pairs_R2(:,2)-all_pairs_R2(:,1));
testing2 = (all_pairs_R2(:,2)-all_pairs_R2(:,1))./all_pairs_R2(:,1);
a = find(testing1 == max(testing1));
b = all_pairs_counter{a,2};

%%% claim one: where ME2/MSE is high, we see the greatest improvement going
% from CTM to CTM + Obs

%%% claim two: 1) where MSB is low and 2) where we are far away from 
% stations, we see the greatest improvement going from Obs to CTM + Obs

% saving all the results
save(sprintf('matfiles/antiperform_%dkm.mat',floor(forceddist./1000)), ...
    'all_R2','all_msBias','all_counter','all_sumidx','all_pairs_sumidx', ...
    'all_pairs_R2','all_pairs_msBias','all_pairs_counter','idx_all_str');

end