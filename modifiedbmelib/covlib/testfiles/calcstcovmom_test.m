function calcstcovmom_test()

sptlcoord = [0,0;0,1;0,2;0,2.9;0,4;0,5];
tempcoord = [1,2,2.9,4];

sdistmat = calcdistmateuclid(sptlcoord,sptlcoord);
tdistmat = calcdistmattemp(tempcoord,tempcoord);

val = [[1;2;3;2;3;1],[2;2;4;3;4;2],[2;2;3;3;3;1],[-1;-2;7;2;1;8]];

vcoord = ones(size(sptlcoord,1)*size(tempcoord,2),3);
for i = 1:size(tempcoord,2)
    vcoord(size(sptlcoord,1)*(i-1)+1:size(sptlcoord,1)*i,:) = ...
       [sptlcoord,tempcoord(i)*ones(size(sptlcoord,1),1)];
end

vval = val(:);

vcoord = vcoord(~isnan(vval),:);
vval = vval(~isnan(vval));

slagbound = [0,1,2,3,4,5];
tlagbound = [0,1,2,3];

[slag,scovval,snumpair,sweight,tlag,tcovval,tnumpair,tweight] = ...
    calcstcovmom(sdistmat,tdistmat,val,slagbound,tlagbound);

[ds1,dt1,c1,o1]=crosscovarioST2(vcoord,vcoord,vval,vval,slagbound,0);
[ds2,dt2,c2,o2]=crosscovarioST2(vcoord,vcoord,vval,vval,0,tlagbound);

disp('calcstcovmom_test')
disp(['Num Pair (mom)            :',num2str(snumpair)]);
disp(['Num Pair (crosscovarioST2):',num2str([size(val(:),1),o1'])]);
disp(['Num Pair (mom)            :',num2str(tnumpair)]);
disp(['Num Pair (crosscovarioST2):',num2str([size(val(:),1),o2])]);


figure;
subplot(2,1,1)
hold on;
plot(slag,scovval,'-ro');
plot(ds1(:,1),c1(:,1),'--b.');

subplot(2,1,2)
hold on;
plot(tlag,tcovval,'-ro');
plot(dt2(1,:),c2(1,:),'--b.');
