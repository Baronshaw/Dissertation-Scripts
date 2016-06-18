function [] = callseeEffectOfSoft()
% this function will call seeEffectOfSoft.m

for i = 0:100000:1000000
    seeEffectOfSoft(i);
    close all;
end

end