function [] = runall_Cluster_04a()
% this function will run all the data in the 04_mfiles_softdata folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../04_mfiles_softdata

uniyrsOff = [ 1998 1999 ; 1999 1999 ; 1999 2000 ; 2000 2000 ; 2000 2001 ; ...
    2001 2001 ; 2001 2002 ; 2002 2002 ; 2002 2003 ; 2003 2003 ; 2003 2004 ; ...
    2004 2004 ; 2004 2005 ; 2005 2005 ; 2005 2006 ; 2006 2006 ; 2006 2007 ; ...
    2007 2007 ; 2007 2008 ; 2008 2008 ; 2008 2009 ; 2009 2009 ; 2009 2010 ; ...
    2010 2010 ; 2010 2011 ; 2011 2011 ; 2011 2012 ];
len = size(uniyrsOff,1);

D = cell(len,1);
Obsg = cell(len,1);
CTMg = cell(len,1);
cMSObs = cell(len,1);
tMEO = cell(len,1);
cMSCTM = cell(len,1);
tMECTM = cell(len,1);
distCTMv = cell(len,1);
yrmodaCTMv = cell(len,1);
dailyCTMv = cell(len,1);

% distances for all the soft data 
for i = 1:len
    [D{i,1} Obsg{i,1} CTMg{i,1} cMSObs{i,1} tMEO{i,1} cMSCTM{i,1} ...
        tMECTM{i,1} distCTMv{i,1} yrmodaCTMv{i,1} dailyCTMv{i,1}] = prepSoft(uniyrsOff(i,:));
end

save('../matfiles/prepSoft.mat','D','Obsg','CTMg','cMSObs','tMEO','cMSCTM', ...
    'tMECTM','distCTMv','yrmodaCTMv','dailyCTMv');

end