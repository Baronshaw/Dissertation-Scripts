function [] = table_MP_cmaq()
% create table of results created from get_MP_cmaq and get_MP_cmaqII

% load results
load('matfiles/traditional_performance.mat');

% header string
str1 = 'num pairs,mean obs,mean mod,mean bias,norm bias,norm mean bias,';
str2 = 'frac bias,mean err,norm err,norm mean err,frac err,corr,corr squared,';
str3 = 's bias,ms bias,rms bias,nrms bias,mean bias DIV s bias,mean bias squared DIV ms bias,v bias DIV ms bias,beta1,v obs,v mod';
strval = sprintf('%s%s%s',str1,str2,str3);

% overall
n = [ num mObs mMod mBias nBias nmBias fBias mErr nErr nmErr fErr R R2 sBias msBias rmsBias nrmsBias mDsBias m2DmsBias s2DmsBias beta1 vObs vMod ];
outid = fopen('tables/overall.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/overall.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

% by percentile
n = [ numP' mObsP' mModP' mBiasP' nBiasP' nmBiasP' fBiasP' mErrP' nErrP' nmErrP' fErrP' RP' R2P' sBiasP' msBiasP' rmsBiasP' nrmsBiasP' mDsBiasP' m2DmsBiasP' s2DmsBiasP' beta1P' vObsP' vModP' ];
outid = fopen('tables/percentile.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/percentile.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

% by season
n = [ numS' mObsS' mModS' mBiasS' nBiasS' nmBiasS' fBiasS' mErrS' nErrS' nmErrS' fErrS' RS' R2S' sBiasS' msBiasS' rmsBiasS' nrmsBiasS' mDsBiasS' m2DmsBiasS' s2DmsBiasS' beta1S' vObsS' vModS' ];
outid = fopen('tables/season.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/season.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

% by land use
n = [ numL' mObsL' mModL' mBiasL' nBiasL' nmBiasL' fBiasL' mErrL' nErrL' nmErrL' fErrL' RL' R2L' sBiasL' msBiasL' rmsBiasL' nrmsBiasL' mDsBiasL' m2DmsBiasL' s2DmsBiasL' beta1L' vObsL' vModL' ];
outid = fopen('tables/landuse.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/landuse.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

% by rank
n = [ numR' mObsR' mModR' mBiasR' nBiasR' nmBiasR' fBiasR' mErrR' nErrR' nmErrR' fErrR' RR' R2R' sBiasR' msBiasR' rmsBiasR' nrmsBiasR' mDsBiasR' m2DmsBiasR' s2DmsBiasR' beta1R' vObsR' vModR' ];
outid = fopen('tables/rank.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/rank.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

% by EPA area
n = [ numE' mObsE' mModE' mBiasE' nBiasE' nmBiasE' fBiasE' mErrE' nErrE' nmErrE' fErrE' RE' R2E' sBiasE' msBiasE' rmsBiasE' nrmsBiasE' mDsBiasE' m2DmsBiasE' s2DmsBiasE' beta1E' vObsE' vModE' ];
outid = fopen('tables/eparegion.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/eparegion.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

% load results
load('matfiles/traditional_performanceII.mat');

% header string
strval = sprintf('va1,val2,%s%s%s',str1,str2,str3);

% by percentile/season
val1 = repmat([1:10]',1,4);
val2 = repmat(1:4,10,1);
n = [ val1(:) val2(:) numPS(:) mObsPS(:) mModPS(:) mBiasPS(:) nBiasPS(:) nmBiasPS(:) fBiasPS(:) mErrPS(:) nErrPS(:) nmErrPS(:) fErrPS(:) RPS(:) R2PS(:) sBiasPS(:) msBiasPS(:) rmsBiasPS(:) nrmsBiasPS(:) mDsBiasPS(:) m2DmsBiasPS(:) s2DmsBiasPS(:) beta1PS(:) vObsPS(:) vModPS(:) ];
outid = fopen('tables/percentile_season.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/percentile_season.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

% by percentile/network
val1 = repmat([1:10]',1,4);
val2 = repmat(1:4,10,1);
n = [ val1(:) val2(:) numPN(:) mObsPN(:) mModPN(:) mBiasPN(:) nBiasPN(:) nmBiasPN(:) fBiasPN(:) mErrPN(:) nErrPN(:) nmErrPN(:) fErrPN(:) RPN(:) R2PN(:) sBiasPN(:) msBiasPN(:) rmsBiasPN(:) nrmsBiasPN(:) mDsBiasPN(:) m2DmsBiasPN(:) s2DmsBiasPN(:) beta1PN(:) vObsPN(:) vModPN(:) ];
outid = fopen('tables/percentile_network.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/percentile_network.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

% by season/network
val1 = repmat([1:4]',1,4);
val2 = repmat(1:4,4,1);
n = [ val1(:) val2(:) numSN(:) mObsSN(:) mModSN(:) mBiasSN(:) nBiasSN(:) nmBiasSN(:) fBiasSN(:) mErrSN(:) nErrSN(:) nmErrSN(:) fErrSN(:) RSN(:) R2SN(:) sBiasSN(:) msBiasSN(:) rmsBiasSN(:) nrmsBiasSN(:) mDsBiasSN(:) m2DmsBiasSN(:) s2DmsBiasSN(:) beta1SN(:) vObsSN(:) vModSN(:) ];
outid = fopen('tables/season_network.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/season_network.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

end