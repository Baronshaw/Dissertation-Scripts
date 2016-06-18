function calcdistmateuclid_test()

sptlcoord = [0,0;1,0;1,1;0,1;-1,1;-1,0;-1,-1;0,-1;1,-1];
distmat = calcdistmateuclid(sptlcoord,sptlcoord);

figure;
hold on;
plot(sptlcoord(:,1),sptlcoord(:,2),'ko','MarkerFaceColor','r');
axis equal;
xlim([-1.2,1.2]);
ylim([-1.2,1.2]);

for i = 1:size(sptlcoord,1)
    for j = 1:size(sptlcoord,1)
        if i~=j
            plot([sptlcoord(i,1),sptlcoord(j,1)],...
                 [sptlcoord(i,2),sptlcoord(j,2)],'-r');
            text((sptlcoord(i,1) + sptlcoord(j,1)) / 2,...
                 (sptlcoord(i,2) + sptlcoord(j,2)) / 2,...
                 num2str(distmat(i,j)));
        end
    end
end

sptlcoord1 = [0,0;0,4];
sptlcoord2 = [0,1;0,2;0,3];
distmat = calcdistmateuclid(sptlcoord1,sptlcoord2);
disp(['(x,y)=(0,0)-(0,1):',num2str(distmat(1,1))]);
disp(['(x,y)=(0,0)-(0,2):',num2str(distmat(1,2))]);
disp(['(x,y)=(0,0)-(0,3):',num2str(distmat(1,3))]);
disp(['(x,y)=(0,4)-(0,1):',num2str(distmat(2,1))]);
disp(['(x,y)=(0,4)-(0,2):',num2str(distmat(2,2))]);
disp(['(x,y)=(0,4)-(0,3):',num2str(distmat(2,3))]);
