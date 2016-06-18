function [] = plotCTMLambdas()
% this function will call CTM and Lambdas for every day/month/year

daystoshow = [datenum(2001,1,1):datenum(2002,12,31) datenum(2005,1,1):datenum(2007,12,31)];
daystoshow = [datenum(2001,1,1):datenum(2001,12,31)];
daystoshow = datevec(daystoshow);
daystoshow = daystoshow(:,1:3);
daystoshowdisp = daystoshow(:,1).*10^4 + daystoshow(:,2).*10^2 + daystoshow(:,3);

% plot individual days
for i = 1:length(daystoshow)
    disp(daystoshow(i,:));
    plotCTM(daystoshow(i,:));
    %plotLambdas(daystoshow(i,:),daystoshowdisp(i));
    %plotCTMbias(daystoshow(i,:),daystoshowdisp(i));
    if mod(i,10)==0, close all; end
end

% % plot average months
% daystoshowmonth = unique(daystoshow(:,1:2),'rows');
% daystoshowmonthdisp = daystoshowmonth(:,1).*10^4 + daystoshowmonth(:,2).*10^2 + 1;
% for i = 1:length(daystoshowmonth)
%     disp(daystoshowmonth(i,:));
%     plotCTMavg(daystoshowmonth(i,:),'mon');
%     cd ../05_mfiles_crossvalidation
%     plotLambdasavg(daystoshowmonth(i,:),daystoshowmonthdisp(i),'mon');
%     cd ../05_mfiles_crossvalidation
%     plotCTMbiasavg(daystoshowmonth(i,:),daystoshowmonthdisp(i),'mon');
%     cd ../04_mfiles_softdata
%     if mod(i,10)==0, close all; end
% end
% 
% % plot average years
% daystoshowyear = unique(daystoshow(:,1),'rows');
% daystoshowyeardisp = daystoshowyear(:,1).*10^4 + 1.*10^2 + 1;
% for i = 1:length(daystoshowyear)
%     disp(daystoshowyear(i,:));
%     plotCTMavg(daystoshowyear(i,:),'yr');
%     cd ../05_mfiles_crossvalidation
%     plotLambdasavg(daystoshowyear(i,:),daystoshowyeardisp(i),'yr');
%     cd ../05_mfiles_crossvalidation
%     plotCTMbiasavg(daystoshowyear(i,:),daystoshowyeardisp(i),'yr');
%     cd ../04_mfiles_softdata
%     if mod(i,10)==0, close all; end
% end

end