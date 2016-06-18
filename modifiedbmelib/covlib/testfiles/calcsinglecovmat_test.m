function calcsinglecovmat_test()

distmat = 0:0.01:15;

manual1 = exponentialC(distmat,[0.7,4]);
output1 = calcsinglecovmat(distmat,'exponentialC',[0.7,4]);
if ~isequal(manual1,output1)
    error('calcsinglecovmat error: exponential nugget');
end

figure;
hold on;
plot(distmat,manual1,'-r','LineWidth',2);
plot(distmat,output1,'-.g','LineWidth',2);
