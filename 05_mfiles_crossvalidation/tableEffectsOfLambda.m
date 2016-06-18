function [] = tableEffectsOfLambda()
% this function will create tables of the results of the lambda1 analysis

soft_years = [2001:2002 2005 2006:2007];
for i = 1:length(soft_years)

    % load results
    load(sprintf('help_lambda_%d.mat',soft_years(i)));
    actuallambda_help2 = ~actuallambda_help;
    test1 = nansum(lambda_help2,2);
    
    % header string
    strval = 'lambda1 help,soft help,probability';

    % table for all data daily
    n1 = [1;1;1;1;0;0;0;0];
    n2 = [0;1;2;3;0;1;2;3];
    n3 = [sum(actuallambda_help(test1==0))/length(actuallambda_help) ; ...
        sum(actuallambda_help(test1==1))/length(actuallambda_help) ; ...
        sum(actuallambda_help(test1==2))/length(actuallambda_help) ; ...
        sum(actuallambda_help(test1==3))/length(actuallambda_help) ; ...
        sum(actuallambda_help2(test1==0))/length(actuallambda_help2) ; ...
        sum(actuallambda_help2(test1==1))/length(actuallambda_help2) ; ...
        sum(actuallambda_help2(test1==2))/length(actuallambda_help2) ; ...
        sum(actuallambda_help2(test1==3))/length(actuallambda_help2) ];
    n = [n1 n2 n3];
    outid = fopen(sprintf('prob_lambda_%d.csv',soft_years(i)),'w+');
    fprintf(outid,'%s',strval);
    fclose(outid);
    dlmwrite(sprintf('prob_lambda_%d.csv',soft_years(i)), ...
        n,'delimiter',',','precision',6,'-append','roffset',1)

end

end