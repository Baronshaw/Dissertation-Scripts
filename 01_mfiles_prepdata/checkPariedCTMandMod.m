function [] = checkPariedCTMandMod()
% this function is to check if the modeled data match the paired modeled
% and obs

CTMyears = [2001:2002 2005 2006 2006 2007];
filenames = { 'J3a_b313.36km.std.2001year_PM.conc' ; ...
    'CCTM_v45_cb4_ae3_aq_ebi_yamo_mpi.Base02b.ACONC.2002yearly.combinePM2p5' ; ...
    'combine.2005cs_05b.metv48.camxv53a.36US1_14.2005yearly.ncf' ; 
    'CCTM47_cb05cl_ae5_aq_64pgf52.aac_spr06.12k.ACONC.2006076-118.ncf.combinePM2p5' ; ...
    'CCTM47_cb05cl_sum06.12k.ACONC.2006181-206.ncf.combinePM2p5' ; ...
    'F1.2007yearly.PM2p5.D36.ncf' };

for i = 3:length(CTMyears)
    
    years = CTMyears(i);
    % load paired modeled and obs
%     if years == 2006 & i == 4
%         fid = fopen(sprintf('../datafiles/Paired_%s/%d/sitecmp_%d_12k_spr06_%s.txt', ...
%             'PM2p5',years,years,'PM2p5'));
%     elseif years == 2006 & i == 5
%         fid = fopen(sprintf('../datafiles/Paired_%s/%d/sitecmp_%d_12k_sum06_%s.txt', ...
%             'PM2p5',years,years,'PM2p5'));
%     else
%         fid = fopen(sprintf('../datafiles/Paired_%s/%d/sitecmp_%d_36k_%s.txt', ...
%             'PM2p5',years,years,'PM2p5'));
%     end
%     ObsVMod = textscan(fid,'%f%f%f%f%f%s%s%f%f', 'delimiter', '\t');
%     fclose(fid);
%     SiteID = ObsVMod{1};
%     LatitudePaired = ObsVMod{2};
%     LongitudePaired = ObsVMod{3};
%     ColumnPaired = ObsVMod{4};
%     RowPaired = ObsVMod{5};
%     TimeOnPaired = ObsVMod{6};
%     TimeOffPaired = ObsVMod{7};
%     Obs = ObsVMod{8};
%     Mod = ObsVMod{9};
%     temp = datenum(datevec(TimeOnPaired));
%     TimeMod = temp-min(temp)+1;

    % load CTM
    cd ../datafiles/Modeled_PM2p5
    if years == 2005
        % I must use the times from a different years, because there are
        % problems with 2005
        tflag = ncread('F1.2007yearly.PM2p5.D36.ncf','TFLAG');
    else
        tflag = ncread(filenames{i},'TFLAG');
    end
    if years == 2005,
        CTM = ncread(filenames{i},'PM25_TOT');
    elseif years == 2007
        CTM = ncread(filenames{i},'PM2.5');
    else
        CTM = ncread(filenames{i},'PM25');
    end
    cd ../../01_mfiles_prepdata
    CTM = squeeze(CTM);
    daysCTM = str2num(num2str(tflag(1,1,:) - floor(tflag(1,1,:)./1000).*1000));
    yearCTM = str2num(num2str(floor(tflag(1,1,:)./1000)));
    hoursCTM = str2num(num2str(floor(tflag(2,1,:)./10000)));

    % average all values to daily values
    [alldays uniidx] = unique(daysCTM);
    yrmodaCTM = datevec(datenum(years,0,alldays));
    [r c h] = size(CTM);
    dailyCTMg = NaN*ones(c,r,length(alldays)); % r and c reverse on purpose
    dailyCTMg_test = NaN*ones(c,r,length(alldays)); % r and c reverse on purpose
    for j = 1:length(alldays)
        idx = alldays(j) == daysCTM;
        dailyCTMg(:,:,j) = mean(CTM(:,:,idx),3)';
        
        a_test = find(idx == 1); a_test = a_test + 4;
        [a b len] = size(CTM);
        a_test(a_test>len) = len;
        dailyCTMg_test(:,:,j) = mean(CTM(:,:,a_test),3)';
    end

    % loop through each obs value and pick the corresponding CTM value 
    Modfile = NaN*ones(length(Obs),1);
    Modfile_test = NaN*ones(length(Obs),1);
    for j = 1:length(Obs)  
        Modfile(j) = dailyCTMg(RowPaired(j),ColumnPaired(j),TimeMod(j));  
        Modfile_test(j) = dailyCTMg_test(RowPaired(j),ColumnPaired(j),TimeMod(j));
    end
    
    % compare the two
    figure; hold on;
    plot([0 max(Mod)],[0 max(Mod)],'r-');
    plot(Modfile,Mod,'b.');
    title(sprintf('Mod file and Mod Paired %d',years));
    xlabel('Modfile');
    ylabel('Mod Paired');
    print(gcf,'-painters','-dpng','-r600',sprintf('ModPaired_Modfile_%d.png',i));
    
    figure; hold on;
    plot([0 max(Mod)],[0 max(Mod)],'r-');
    plot(Modfile_test,Mod,'b.');
    title(sprintf('Mod file +4 and Mod Paired %d',years));
    xlabel('Modfile +4');
    ylabel('Mod Paired');
    print(gcf,'-painters','-dpng','-r600',sprintf('ModPaired_Modfile4_%d.png',i));

end

end