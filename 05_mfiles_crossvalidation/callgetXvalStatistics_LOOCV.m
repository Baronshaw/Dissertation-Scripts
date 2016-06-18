function [] = callgetXvalStatistics_LOOCV()

forceddist = 0:100000:1000000;

for i = 0:1
    for j = 1:length(forceddist)
        getXvalStatistics_LOOCV(i,0,1,forceddist(j));
    end
end

end