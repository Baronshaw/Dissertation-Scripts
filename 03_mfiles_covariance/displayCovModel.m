function [] = displayCovModel()
% this function will display the results of the covariance modeling in a file

% load information for each offset/model
meanNow = {'short','intermediate','long','very_long'}; 

%%% in space
% loading all the model names and their parameters
load('../matfiles/covmod_modinfo_joint.mat');

n = 1;
for i = 1:length(meanNow)
    for j = 1:length(modstr)
        
        offstrr{n,1} = meanNow{i};
        modstrr{n,1} = modstr{i,j};
        load(sprintf('../matfiles/covmod_r_%s_%s_joint.mat',meanNow{i},modstr{i,j}));
        r2r{n,1} = gof.rsquare;
        AICr{n,1} = aicval;
        paramr{n,1} = coeffnames(f)';
        paramvalr{n,1} = coeffvalues(f)';
        n = n + 1;
        
    end
end

%%% in time
% loading all the model names and their parameters
load('../matfiles/covmod_modinfo_joint.mat');

n = 1;
for i = 1:length(meanNow)
    for j = 1:length(modstr)
        
        offstrt{n,1} = meanNow{i};
        modstrt{n,1} = modstr{i,j};
        load(sprintf('../matfiles/covmod_t_%s_%s_joint.mat',meanNow{i},modstr{i,j}));
        r2t{n,1} = gof.rsquare;
        AICt{n,1} = aicval;
        paramt{n,1} = coeffnames(f)';
        paramvalt{n,1} = coeffvalues(f)';
        n = n + 1;
        
    end
end

% writing to a .csv file
charparam = cellfun(@anon1,[paramr;paramt],'UniformOutput',false);
charparamval = cellfun(@anon2,[paramvalr;paramvalt],'UniformOutput',false);
offstr = [offstrr;offstrt];
modstr = [modstrr;modstrt];
r2 = [r2r;r2t];
AIC = [AICr;AICt];

fid = fopen('covmodjoint_results.csv','wt');
for i = 1:size(charparam,1)
    fprintf(fid, '%s,%s,%0.2f,%0.2f,%s,%s\n',offstr{i},modstr{i},...
        r2{i},AIC{i},charparam{i},charparamval{i});
end
fclose(fid);

% writing all the AIC values to a .csv file
load('../matfiles/covmod_modinfo_joint.mat');
fid = fopen('covmodjoint_results_AIC.csv','wt');
fprintf(fid,'space');
for i = 1:length(meanNow)
    fprintf(fid,',%s',meanNow{i});
end
fprintf(fid,'\n');
for i = 1:length(modstr)
    fprintf(fid,'%s,',modstr{1,i});
    for j = 1:length(meanNow)
        load(sprintf('../matfiles/covmod_r_%s_%s_joint.mat',meanNow{j},modstr{j,i}));
        fprintf(fid,'%0.2f,',aicval);
    end
    fprintf(fid,'\n');
end
fprintf(fid,'time');
for i = 1:length(meanNow)
    fprintf(fid,',%s',meanNow{i});
end
fprintf(fid,'\n');
for i = 1:length(modstr)
    fprintf(fid,'%s,',modstr{1,i});
    for j = 1:length(meanNow)
        load(sprintf('../matfiles/covmod_t_%s_%s_joint.mat',meanNow{j},modstr{j,i}));
        fprintf(fid,'%0.2f,',aicval);
    end
    fprintf(fid,'\n');
end
fclose(fid);

% writing all the r2 values to a .csv file
load('../matfiles/covmod_modinfo_joint.mat');
fid = fopen('covmodjoint_results_r2.csv','wt');
fprintf(fid,'space');
for i = 1:length(meanNow)
    fprintf(fid,',%s',meanNow{i});
end
fprintf(fid,'\n');
for i = 1:length(modstr)
    fprintf(fid,'%s,',modstr{1,i});
    for j = 1:length(meanNow)
        load(sprintf('../matfiles/covmod_r_%s_%s_joint.mat',meanNow{j},modstr{j,i}));
        fprintf(fid,'%0.4f,',gof.rsquare);
    end
    fprintf(fid,'\n');
end
fprintf(fid,'time');
for i = 1:length(meanNow)
    fprintf(fid,',%s',meanNow{i});
end
fprintf(fid,'\n');
for i = 1:length(modstr)
    fprintf(fid,'%s,',modstr{1,i});
    for j = 1:length(meanNow)
        load(sprintf('../matfiles/covmod_t_%s_%s_joint.mat',meanNow{j},modstr{j,i}));
        fprintf(fid,'%0.4f,',gof.rsquare);
    end
    fprintf(fid,'\n');
end
fclose(fid);

end

function [stringA] = anon1(blah)
    stringA = char(blah{1});
    for n = 2:length(blah)
        stringA = [stringA sprintf('/%s',blah{n})];
    end

end

function [stringA] = anon2(blah)
    stringA = sprintf('%0.2f',blah(1));
    for n = 2:length(blah)
        stringA = [stringA sprintf('/%0.2f',blah(n))];
    end
end
