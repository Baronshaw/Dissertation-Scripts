function [D Obsg CTMg cMSObs tMEO cMSCTM tMECTM distCTMv yrmodaCTMv dailyCTMv] = prepSoft(uniyrsOff)
% this function will prep the soft data by finding the distances between
% the CTM data and Observational data


if exist(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1))) | ...
        exist(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(2)))
    
    if exist(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1))) & ...
            ~exist(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(2)))
        % only first exists

        load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1)));
        load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(1)))

    elseif ~exist(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1))) & ...
            exist(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(2)))
        % only second exists

        load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(2)));
        load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(2)))

    elseif exist(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1))) & ...
            exist(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(2)))
        % both exist

        if uniyrsOff(1) == uniyrsOff(2) % both years equal
            if uniyrsOff(1) == 2006 % exception for 2006
                
                load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1)));
                load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(1)))
                coordObs1 = coordObs; Obs1 = Obs; Mod1 = Mod; yrmodaObs1 = yrmodaObs;
                distCTMv1=distCTMv; yrmodaCTMv1=yrmodaCTMv; dailyCTMv1=dailyCTMv;
                load(sprintf('../matfiles/prepCTMandObs_%dsum.mat',uniyrsOff(2)));
                load(sprintf('../matfiles/prepCTM_%dsum.mat',uniyrsOff(2)));
                coordObs2 = coordObs; Obs2 = Obs; Mod2 = Mod; yrmodaObs2 = yrmodaObs;
                distCTMv2=distCTMv; yrmodaCTMv2=yrmodaCTMv; dailyCTMv2=dailyCTMv;
                coordObs=[coordObs1;coordObs2]; Obs=[Obs1;Obs2]; Mod=[Mod1;Mod2]; yrmodaObs=[yrmodaObs1;yrmodaObs2];
                distCTMv=[distCTMv1;distCTMv2]; yrmodaCTMv=[yrmodaCTMv1;yrmodaCTMv2]; dailyCTMv=[dailyCTMv1;dailyCTMv2];
                
            else
                load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1)));
                load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(1)))
            end
        else % both years different
            
            if uniyrsOff(1) == 2006
                
                load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1)));
                load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(1)))
                coordObs1 = coordObs; Obs1 = Obs; Mod1 = Mod; yrmodaObs1 = yrmodaObs;
                distCTMv1=distCTMv; yrmodaCTMv1=yrmodaCTMv; dailyCTMv1=dailyCTMv;
                load(sprintf('../matfiles/prepCTMandObs_%dsum.mat',uniyrsOff(1)));
                load(sprintf('../matfiles/prepCTM_%dsum.mat',uniyrsOff(1)));
                coordObs2 = coordObs; Obs2 = Obs; Mod2 = Mod; yrmodaObs2 = yrmodaObs;
                distCTMv2=distCTMv; yrmodaCTMv2=yrmodaCTMv; dailyCTMv2=dailyCTMv;
                load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(2)));
                load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(2)));
                coordObs3 = coordObs; Obs3 = Obs; Mod3 = Mod; yrmodaObs3 = yrmodaObs;
                distCTMv3=distCTMv; yrmodaCTMv3=yrmodaCTMv; dailyCTMv3=dailyCTMv;
                coordObs=[coordObs1;coordObs2;coordObs3]; Obs=[Obs1;Obs2;Obs3]; Mod=[Mod1;Mod2;Mod3]; yrmodaObs=[yrmodaObs1;yrmodaObs2;yrmodaObs3];
                distCTMv=[distCTMv1;distCTMv2;distCTMv3]; yrmodaCTMv=[yrmodaCTMv1;yrmodaCTMv2;yrmodaCTMv3]; dailyCTMv=[dailyCTMv1;dailyCTMv2;dailyCTMv3];

                
            elseif uniyrsOff(2) == 6
                
                load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1)));
                load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(1)));
                coordObs1 = coordObs; Obs1 = Obs; Mod1 = Mod; yrmodaObs1 = yrmodaObs;
                distCTMv1=distCTMv; yrmodaCTMv1=yrmodaCTMv; dailyCTMv1=dailyCTMv;
                load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(2)));
                load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(2)))
                coordObs2 = coordObs; Obs2 = Obs; Mod2 = Mod; yrmodaObs2 = yrmodaObs;
                distCTMv2=distCTMv; yrmodaCTMv2=yrmodaCTMv; dailyCTMv2=dailyCTMv;
                load(sprintf('../matfiles/prepCTMandObs_%dsum.mat',uniyrsOff(2)));
                load(sprintf('../matfiles/prepCTM_%dsum.mat',uniyrsOff(2)));
                coordObs3 = coordObs; Obs3 = Obs; Mod3 = Mod; yrmodaObs3 = yrmodaObs;
                distCTMv3=distCTMv; yrmodaCTMv3=yrmodaCTMv; dailyCTMv3=dailyCTMv;                
                coordObs=[coordObs1;coordObs2;coordObs3]; Obs=[Obs1;Obs2;Obs3]; Mod=[Mod1;Mod2;Mod3]; yrmodaObs=[yrmodaObs1;yrmodaObs2;yrmodaObs3];
                distCTMv=[distCTMv1;distCTMv2;distCTMv3]; yrmodaCTMv=[yrmodaCTMv1;yrmodaCTMv2;yrmodaCTMv3]; dailyCTMv=[dailyCTMv1;dailyCTMv2;dailyCTMv3];
                
            else

                load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(1)));
                load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(1)))
                coordObs1 = coordObs; Obs1 = Obs; Mod1 = Mod; yrmodaObs1 = yrmodaObs;
                distCTMv1=distCTMv; yrmodaCTMv1=yrmodaCTMv; dailyCTMv1=dailyCTMv;
                load(sprintf('../matfiles/prepCTMandObs_%d.mat',uniyrsOff(2)));
                load(sprintf('../matfiles/prepCTM_%d.mat',uniyrsOff(2)));
                coordObs2 = coordObs; Obs2 = Obs; Mod2 = Mod; yrmodaObs2 = yrmodaObs;
                distCTMv2=distCTMv; yrmodaCTMv2=yrmodaCTMv; dailyCTMv2=dailyCTMv;
                coordObs=[coordObs1;coordObs2]; Obs=[Obs1;Obs2]; Mod=[Mod1;Mod2]; yrmodaObs=[yrmodaObs1;yrmodaObs2];
                distCTMv=[distCTMv1;distCTMv2]; yrmodaCTMv=[yrmodaCTMv1;yrmodaCTMv2]; dailyCTMv=[dailyCTMv1;dailyCTMv2];

            end
            
        end

    end

        [Obsg,cMSObs,tMEO,nanratio]=valstv2stg_new([coordObs(:,1) coordObs(:,2) yrmodaObs],Obs);
        [CTMg,dummy,dummy,nanratio]=valstv2stg_new([coordObs(:,1) coordObs(:,2) yrmodaObs],Mod);
        yrmodaCTMv = yrmodaCTMv(:,1)*10000 + yrmodaCTMv(:,2)*100 + yrmodaCTMv(:,3);
        [dummy,cMSCTM,tMECTM,nanratio]=valstv2stg_new([distCTMv yrmodaCTMv],dailyCTMv);

        Y = [cMSObs(:,1) cMSObs(:,2)]; 
        Y = Y';
        X = [cMSCTM(:,1) cMSCTM(:,2)]';
        D = bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y);
   
else % neither files exists
    
    D=[]; Obsg=[]; CTMg=[]; cMSObs=[]; tMEO=[]; cMSCTM=[]; tMECTM=[];
    distCTMv=[]; yrmodaCTMv=[]; dailyCTMv=[];
   
end

end