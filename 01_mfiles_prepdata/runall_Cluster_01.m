function [] = runall_Cluster_01()
% this function will run all the data in the 01_mfiles_prepdata folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../01_mfiles_prepdata

% files names and years of interest
years = 1999:2011;
CTMyears = [2001:2002 2005 2006 2006 2007];
CTMyears = [2001:2002];
CTMyears = 2005;
deltaT = 365;
filenames = { 'J3a_b313.36km.std.2001year_PM.conc' ; ...
    'CCTM_v45_cb4_ae3_aq_ebi_yamo_mpi.Base02b.ACONC.2002yearly.combinePM2p5' ; ...
    'combine.2005cs_05b.metv48.camxv53a.36US1_14.2005yearly.ncf' ; 
    'CCTM47_cb05cl_ae5_aq_64pgf52.aac_spr06.12k.ACONC.2006076-118.ncf.combinePM2p5' ; ...
    'CCTM47_cb05cl_sum06.12k.ACONC.2006181-206.ncf.combinePM2p5' ; ...
    'F1.2007yearly.PM2p5.D36.ncf' };
filenames = { 'combine.2005cs_05b.metv48.camxv53a.36US1_14.2005yearly.ncf' ; 
    'CCTM47_cb05cl_ae5_aq_64pgf52.aac_spr06.12k.ACONC.2006076-118.ncf.combinePM2p5' ; ...
    'CCTM47_cb05cl_sum06.12k.ACONC.2006181-206.ncf.combinePM2p5' ; ...
    'F1.2007yearly.PM2p5.D36.ncf' };

% prep CTM data
for i = 1:length(CTMyears) 
    tic
    disp(CTMyears(i));
    [whiproj nad83 distCTM coordCTM yrmodaCTM dailyCTMg distCTMv ...
    yrmodaCTMv dailyCTMv] = prepCTM(filenames{i},CTMyears(i));
    if CTMyears(i) == 2006 && i == 5
        save(sprintf('../matfiles/prepCTM_%dsum.mat',CTMyears(i)),'whiproj','nad83','distCTM',...
            'coordCTM','yrmodaCTM','dailyCTMg','distCTMv','yrmodaCTMv','dailyCTMv'); 
    else
        save(sprintf('../matfiles/prepCTM_%d.mat',CTMyears(i)),'whiproj','nad83','distCTM',...
            'coordCTM','yrmodaCTM','dailyCTMg','distCTMv','yrmodaCTMv','dailyCTMv'); % added 9/14/2014
    end
    toc
end

% prep Obs data
for i = 1:length(years) 
    tic
    disp(years(i));
    [coordObs Obs Mod yrmodaObs] = prepObs(years(i));
    save(sprintf('../matfiles/prepObs_%d.mat',years(i)),'coordObs','Obs','Mod','yrmodaObs');  
    toc
end

% prep Mod and Obs data
for i = 1:length(CTMyears)
    tic
    disp(CTMyears(i));
    [coordObs Obs Mod yrmodaObs] = prepCTMandObs(CTMyears(i));
    if CTMyears(i) == 2006 && i == 5
        save(sprintf('../matfiles/prepCTMandObs_%dsum.mat',CTMyears(i)),'coordObs','Obs','Mod','yrmodaObs');
    else
        save(sprintf('../matfiles/prepCTMandObs_%d.mat',CTMyears(i)),'coordObs','Obs','Mod','yrmodaObs');  
    end
    toc
end

end