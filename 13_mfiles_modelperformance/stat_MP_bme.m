function [] = stat_MP_bme(forceddist,soft,constant,gauss,constr,gaussstr)
% this function will calculate the cross validation statistics for the 
% LOOCV cross validation

% Traditional Measures of Model Performance Include:
% from http://www.epa.gov/cleanairactbenefits/jan10/812_AQM_mpe_memo.pdf
% 1) number of paired modeled and obs
% 2) mean obs value
% 3) mean modeled value
% 4) mean bias
% 5) normalized bias
% 6) normalized mean bias
% 7) fractional bias
% 8) mean error
% 9) normalized error
% 10) normalized mean error
% 11) fractional error 
% 12) correlation
% proposed new measures
% 13) correlation squared
% 14) standard bias
% 15) mean squared bias
% 16) root mean squared bias
% 17) normalized root mean squared bias
% 18) mean bias/standard bias
% 19) mean bias squared/mean squared bias
% 20) variance of bias/mean squared bias
% other measures
% 21) beta1
% 22) variance of obs
% 23) variance of mod

% In BME language
% 1) num of paired modeled and obs
% 2) mean obs 
% 3) mean mod 
% 4) mean error 
% 5) mean normalized error
% 6) normalized mean error 
% 7) fractional error 
% 8) mean absolute error
% 9) mean normalized absolute error 
% 10) normalized mean absolute error 
% 11) fractional absolute error
% 12) correlation 
% 13) correlation squared 
% 14) standard error 
% 15) mean squared error
% 16) root mean squared error 
% 17) normalized root mean squared error 
% 18) mean error/standard error 
% 19) mean error squared/mean squared error
% 20) variance of errors/mean squared error 
% 21) beta1 
% 22) var obs 
% 23) var mod

% parameters
if nargin < 1, forceddist = 0; end
if nargin < 2, soft = 1; end 
if nargin < 3, constant = 0; end
if nargin < 4, gauss = 1; end
if nargin < 5, constr = '_long'; end
if nargin < 6, gausstr = '_gauss'; end

% parameters
constant = 0; gauss = 1; constr = '_long'; gaussstr = '_gauss'; 

for i = 0:2

    soft = i;
    if soft == 0, softstr = '_nosoft'; 
    elseif soft == 1, softstr = '_soft'; 
    end
    
    if i == 0 | i == 1
        disp('here1');
        % loading results
        if forceddist == 0
            load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone__20150408.mat', ...
                softstr,constr,gaussstr,floor(forceddist./1000)));
        else
            load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone_20150408.mat', ...
                softstr,constr,gaussstr,floor(forceddist./1000)));
        end
        zk = cell2mat(zk_madd);
        zh = cell2mat(zh_Xval);
        vk = cell2mat(vk);
        ck = cell2mat(ckXval);

        % just for 2001
        datez = datevec(ck(:,3));
        idx = datez(:,1) == 2001 & ~isnan(zk);
        zk = zk(idx); zh = zh(idx); ck = ck(idx,:); vk = vk(idx);

        % overall
        Obs = zh; Mod = zk;
        num = size(Obs,1); 
        mObs = mean(Obs); 
        mMod = mean(Mod); 
        mBias = mean(Mod-Obs);
        idxff = Obs > 0; % adding a little fudge factor
        nBias = 100.*(1./num).*sum((Mod(idxff)-Obs(idxff))./Obs(idxff)); 
        nmBias = sum(Mod-Obs)./sum(Obs);
        fBias = 100.*(1./num).*sum((Mod-Obs)./(0.5.*(Mod+Obs))); 
        mErr = mean(abs(Mod-Obs)); 
        nErr = 100.*(1./num).*sum(abs(Mod(idxff)-Obs(idxff))./Obs(idxff));
        nmErr = sum(abs(Mod-Obs))./sum(Obs); 
        fErr = 100.*(1./num).*sum(abs(Mod-Obs)./(0.5.*(Mod+Obs))); 
        R = corr(Mod,Obs); 
        R2 = R.^2;
        sBias = std(Mod-Obs);
        msBias = mean((Mod-Obs).^2);
        rmsBias = sqrt(mean((Mod-Obs).^2));
        nrmsBias = sqrt(mean((Mod-Obs).^2))./mean(Obs);
        mDsBias = mBias./sBias;
        m2DmsBias = (mean(Mod-Obs)).^2./mean((Mod-Obs).^2);
        s2DmsBias = var(Mod-Obs)./mean((Mod-Obs).^2);
        dummy = polyfit(Obs,Mod,1);
        beta1 = dummy(1);
        vObs = var(Obs);
        vMod = var(Mod);

        % by station
        uniPMM = unique(ck(:,1:2),'rows');
        for j = 1:length(uniPMM)
            idx = ck(:,1) == uniPMM(j,1) & ck(:,2) == uniPMM(j,2);
            numM(j) = size(Obs(idx),1); 
            mObsM(j) = mean(Obs(idx)); 
            mModM(j) = mean(Mod(idx)); 
            mBiasM(j) = mean(Mod(idx)-Obs(idx));
            idxff = Obs > 0; % adding a little fudge factor
            nBiasM(j) = 100.*(1./numM(j)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
            nmBiasM(j) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
            fBiasM(j) = 100.*(1./numM(j)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
            mErrM(j) = mean(abs(Mod(idx)-Obs(idx))); 
            nErrM(j) = 100.*(1./numM(j)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
            nmErrM(j) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
            fErrM(j) = 100.*(1./numM(j)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
            RM(j) = corr(Mod(idx),Obs(idx));  
            R2M(j) = RM(j).^2;
            sBiasM(j) = std(Mod(idx)-Obs(idx));
            msBiasM(j) = mean((Mod(idx)-Obs(idx)).^2);
            rmsBiasM(j) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
            nrmsBiasM(j) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
            mDsBiasM(j) = mBiasM(j)./sBiasM(j);
            m2DmsBiasM(j) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
            s2DmsBiasM(j) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
            dummy = polyfit(Obs(idx),Mod(idx),1);
            beta1M(j) = dummy(1);
            vObsM(j) = var(Obs(idx));
            vModM(j) = var(Mod(idx));
        end

        % save results
        save(sprintf('matfiles/traditional_performance_bme%s_%dkm.mat',softstr,floor(forceddist./1000)), ...
            'num','mObs','mMod','mBias','nBias','nmBias','fBias','mErr','nErr','nmErr','fErr','R', ... 
            'R2','sBias','msBias','rmsBias','nrmsBias','mDsBias','m2DmsBias','s2DmsBias','beta1','vObs','vMod', ...
            'uniPMM','numM','mObsM','mModM','mBiasM','nBiasM','nmBiasM','fBiasM','mErrM','nErrM','nmErrM','fErrM','RM', ... 
            'R2M','sBiasM','msBiasM','rmsBiasM','nrmsBiasM','mDsBiasM','m2DmsBiasM','s2DmsBiasM','beta1M','vObsM','vModM');
    
    elseif i == 2
        disp('here2');
        load(sprintf('matfiles/traditional_performance_bme%s_%dkm.mat','_nosoft',floor(forceddist./1000)));
        num1 = num; mObs1 = mObs; mMod1 = mMod; mBias1 = mBias; nBias1 = nBias;
        nmBias1 = nmBias; fBias1 = fBias; mErr1 = mErr; nErr1 = nErr; nmErr1 = nmErr;
        fErr1 = fErr; R1 = R; R21 = R2; sBias1 = sBias; msBias1 = msBias; rmsBias1 = rmsBias;
        nrmsBias1 = nrmsBias; mDsBias1 = mDsBias; m2DmsBias1 = m2DmsBias;
        s2DmsBias1 = s2DmsBias; beta11 = beta1; vObs1 = vObs; vMod1 = vMod;
        
        uniPMM1 = uniPMM; numM1 = numM; mObsM1 = mObsM; mModM1 = mModM; mBiasM1 = mBiasM; nBiasM1 = nBiasM;
        nmBiasM1 = nmBiasM; fBiasM1 = fBiasM; mErrM1 = mErrM; nErrM1 = nErrM; nmErrM1 = nmErrM;
        fErrM1 = fErrM; RM1 = RM; R2M1 = R2M; sBiasM1 = sBiasM; msBiasM1 = msBiasM; rmsBiasM1 = rmsBiasM;
        nrmsBiasM1 = nrmsBiasM; mDsBiasM1 = mDsBiasM; m2DmsBiasM1 = m2DmsBiasM;
        s2DmsBiasM1 = s2DmsBiasM; beta1M1 = beta1M; vObsM1 = vObsM; vModM1 = vModM;
        
        load(sprintf('matfiles/traditional_performance_bme%s_%dkm.mat','_soft',floor(forceddist./1000)));
        num2 = num; mObs2 = mObs; mMod2 = mMod; mBias2 = mBias; nBias2 = nBias;
        nmBias2 = nmBias; fBias2 = fBias; mErr2 = mErr; nErr2 = nErr; nmErr2 = nmErr;
        fErr2 = fErr; R2a = R; R22 = R2; sBias2 = sBias; msBias2 = msBias; rmsBias2 = rmsBias;
        nrmsBias2 = nrmsBias; mDsBias2 = mDsBias; m2DmsBias2 = m2DmsBias;
        s2DmsBias2 = s2DmsBias; beta12 = beta1; vObs2 = vObs; vMod2 = vMod;
        
        uniPMM2 = uniPMM; numM2 = numM; mObsM2 = mObsM; mModM2 = mModM; mBiasM2 = mBiasM; nBiasM2 = nBiasM;
        nmBiasM2 = nmBiasM; fBiasM2 = fBiasM; mErrM2 = mErrM; nErrM2 = nErrM; nmErrM2 = nmErrM;
        fErrM2 = fErrM; RM2 = RM; R2M2 = R2M; sBiasM2 = sBiasM; msBiasM2 = msBiasM; rmsBiasM2 = rmsBiasM;
        nrmsBiasM2 = nrmsBiasM; mDsBiasM2 = mDsBiasM; m2DmsBiasM2 = m2DmsBiasM;
        s2DmsBiasM2 = s2DmsBiasM; beta1M2 = beta1M; vObsM2 = vObsM; vModM2 = vModM;
        
        num = 100.*(num2-num1)./num1;
        mObs = 100.*(mObs2-mObs1)./mObs1;
        mMod = 100.*(mMod2-mMod1)./mMod1;
        mBias = 100.*(mBias2-mBias1)./mBias1;
        nBias = 100.*(nBias2-nBias1)./nBias1;
        nmBias = 100.*(nmBias2-nmBias1)./nmBias1;
        fBias = 100.*(fBias2-fBias1)./fBias1;
        mErr = 100.*(mErr2-mErr1)./mErr1;
        nErr = 100.*(nErr2-nErr1)./nErr1;
        nmErr = 100.*(nmErr2-nmErr1)./nmErr1;
        fErr = 100.*(fErr2-fErr1)./fErr1;
        R = 100.*(R2a-R1)./R1;
        R2 = 100.*(R22-R21)./(R21);
        sBias = 100.*(sBias2-sBias1)./sBias1;
        msBias = 100.*(msBias2-msBias1)./msBias1;
        rmsBias = 100.*(rmsBias2-rmsBias1)./rmsBias1;
        nrmsBias = 100.*(nrmsBias2-nrmsBias1)./nrmsBias1;
        mDsBias = 100.*(mDsBias2-mDsBias1)./mDsBias1;
        m2DmsBias = 100.*(m2DmsBias2-m2DmsBias1)./m2DmsBias1;
        s2DmsBias = 100.*(s2DmsBias2-s2DmsBias1)./s2DmsBias1;
        beta1 = 100.*(beta12-beta11)./beta11;
        vObs = 100.*(vObs2-vObs1)./vObs1;
        vMod = 100.*(vMod2-vMod1)./vMod1;
        uniPMM = uniPMM1;
        numM = 100.*(numM2-numM1)./numM1;
        mObsM = 100.*(mObsM2-mObsM1)./mObsM1;
        mModM = 100.*(mModM2-mModM1)./mModM1;
        mBiasM = 100.*(mBiasM2-mBiasM1)./mBiasM1;
        nBiasM = 100.*(nBiasM2-nBiasM1)./nBiasM1;
        nmBiasM = 100.*(nmBiasM2-nmBiasM1)./nmBiasM1;
        fBiasM = 100.*(fBiasM2-fBiasM1)./fBiasM1;
        mErrM = 100.*(mErrM2-mErrM1)./mErrM1;
        nErrM = 100.*(nErrM2-nErrM1)./nErrM1;
        nmErrM = 100.*(nmErrM2-nmErrM1)./nmErrM1;
        fErrM = 100.*(fErrM2-fErrM1)./fErrM1;
        RM = 100.*(RM2-RM1)./RM1;
        R2M = 100.*(R2M2-R2M1)./(R2M1);
        sBiasM = 100.*(sBiasM2-sBiasM1)./sBiasM1;
        msBiasM = 100.*(msBiasM2-msBiasM1)./msBiasM1;
        rmsBiasM = 100.*(rmsBiasM2-rmsBiasM1)./rmsBiasM1;
        nrmsBiasM = 100.*(nrmsBiasM2-nrmsBiasM1)./nrmsBiasM1;
        mDsBiasM = 100.*(mDsBiasM2-mDsBiasM1)./mDsBiasM1;
        m2DmsBiasM = 100.*(m2DmsBiasM2-m2DmsBiasM1)./m2DmsBiasM1;
        s2DmsBiasM = 100.*(s2DmsBiasM2-s2DmsBiasM1)./s2DmsBiasM1;
        beta1M = 100.*(beta1M2-beta1M1)./beta1M1;
        vObsM = 100.*(vObsM2-vObsM1)./vObsM1;
        vModM = 100.*(vModM2-vModM1)./vModM1;
        
        % save results
        save(sprintf('matfiles/traditional_performance_bme_reldiff_%dkm.mat',floor(forceddist./1000)), ...
            'num','mObs','mMod','mBias','nBias','nmBias','fBias','mErr','nErr','nmErr','fErr','R', ... 
            'R2','sBias','msBias','rmsBias','nrmsBias','mDsBias','m2DmsBias','s2DmsBias','beta1','vObs','vMod', ...
            'uniPMM','numM','mObsM','mModM','mBiasM','nBiasM','nmBiasM','fBiasM','mErrM','nErrM','nmErrM','fErrM','RM', ... 
            'R2M','sBiasM','msBiasM','rmsBiasM','nrmsBiasM','mDsBiasM','m2DmsBiasM','s2DmsBiasM','beta1M','vObsM','vModM');
    
        num = (num2-num1);
        mObs = (mObs2-mObs1);
        mMod = (mMod2-mMod1);
        mBias = (mBias2-mBias1);
        nBias = (nBias2-nBias1);
        nmBias = (nmBias2-nmBias1);
        fBias = (fBias2-fBias1);
        mErr = (mErr2-mErr1);
        nErr = (nErr2-nErr1);
        nmErr = (nmErr2-nmErr1);
        fErr = (fErr2-fErr1);
        R = (R2a-R1);
        R2 = (R22-R21);
        sBias = (sBias2-sBias1);
        msBias = (msBias2-msBias1);
        rmsBias = (rmsBias2-rmsBias1);
        nrmsBias = (nrmsBias2-nrmsBias1);
        mDsBias = (mDsBias2-mDsBias1);
        m2DmsBias = (m2DmsBias2-m2DmsBias1);
        s2DmsBias = (s2DmsBias2-s2DmsBias1);
        beta1 = (beta12-beta11);
        vObs = (vObs2-vObs1);
        vMod = (vMod2-vMod1);
        uniPMM = uniPMM1;
        numM = (numM2-numM1);
        mObsM = (mObsM2-mObsM1);
        mModM = (mModM2-mModM1);
        mBiasM = (mBiasM2-mBiasM1);
        nBiasM = (nBiasM2-nBiasM1);
        nmBiasM = (nmBiasM2-nmBiasM1);
        fBiasM = (fBiasM2-fBiasM1);
        mErrM = (mErrM2-mErrM1);
        nErrM = (nErrM2-nErrM1);
        nmErrM = (nmErrM2-nmErrM1);
        fErrM = (fErrM2-fErrM1);
        RM = (RM2-RM1);
        R2M = (R2M2-R2M1);
        sBiasM = (sBiasM2-sBiasM1);
        msBiasM = (msBiasM2-msBiasM1);
        rmsBiasM = (rmsBiasM2-rmsBiasM1);
        nrmsBiasM = (nrmsBiasM2-nrmsBiasM1);
        mDsBiasM = (mDsBiasM2-mDsBiasM1);
        m2DmsBiasM = (m2DmsBiasM2-m2DmsBiasM1);
        s2DmsBiasM = (s2DmsBiasM2-s2DmsBiasM1);
        beta1M = (beta1M2-beta1M1);
        vObsM = (vObsM2-vObsM1);
        vModM = (vModM2-vModM1);
        
        % save results
        save(sprintf('matfiles/traditional_performance_bme_diff_%dkm.mat',floor(forceddist./1000)), ...
            'num','mObs','mMod','mBias','nBias','nmBias','fBias','mErr','nErr','nmErr','fErr','R', ... 
            'R2','sBias','msBias','rmsBias','nrmsBias','mDsBias','m2DmsBias','s2DmsBias','beta1','vObs','vMod', ...
            'uniPMM','numM','mObsM','mModM','mBiasM','nBiasM','nmBiasM','fBiasM','mErrM','nErrM','nmErrM','fErrM','RM', ... 
            'R2M','sBiasM','msBiasM','rmsBiasM','nrmsBiasM','mDsBiasM','m2DmsBiasM','s2DmsBiasM','beta1M','vObsM','vModM');
            
    end
    
end

end