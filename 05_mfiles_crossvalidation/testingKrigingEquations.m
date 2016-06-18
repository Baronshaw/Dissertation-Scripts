function [] = testingKrigingEquations()
% this stupid function will be used to justify my explanations of how sills
% and ranges effect BME variance. I shouldn't have to do this, but I know
% this is the only way that my suggestions will be considered.

% baseline example
ch = [0;0.5;3.1];
zh = [20;14;18];
model = 'exponentialC';
param = [1 1];
nhmax = 3;
dmax = 100;
order = 0;
ck = 1.7;
[zk,vk]=kriging(ck,ch,zh,model,param,nhmax,dmax,order);
% zk = 17.2969; vk = 1.3472;

% decrease variance
param = [0.5 1];
[zk,vk]=kriging(ck,ch,zh,model,param,nhmax,dmax,order);
% zk = 17.2969; vk = 0.6736;

% increase variance
param = [2 1];
[zk,vk]=kriging(ck,ch,zh,model,param,nhmax,dmax,order);
% zk = 17.2969; vk = 2.6944;

% increasing range
param = [1 2];
[zk,vk]=kriging(ck,ch,zh,model,param,nhmax,dmax,order);
% zk = 16.9524; vk = 1.1808;

% decreasing range
param = [1 0.5];
[zk,vk]=kriging(ck,ch,zh,model,param,nhmax,dmax,order);
% zk = 17.3419; vk = 1.3435;

% large range
param = [1 10];
[zk,vk]=kriging(ck,ch,zh,model,param,nhmax,dmax,order);
% zk = 16.0011; vk = 0.3728;

% small range
param = [1 0.1];
[zk,vk]=kriging(ck,ch,zh,model,param,nhmax,dmax,order);
% zk = 17.3333; vk = 1.3333;

end