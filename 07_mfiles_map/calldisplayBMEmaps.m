function [] = calldisplayBMEmaps()

% for i = 0:1
%     for j = 2001:2002
%         for k = 1:12
%             displayBMEmaps([j k 15],i)
%             close all;
%         end
%     end
% end
matlabpool open 12
for i = 0:1
    for j = 1:12
        tic
        BMEmaps([2001 j 1],i,0,1);
        toc
    end
end
matlabpool close

end