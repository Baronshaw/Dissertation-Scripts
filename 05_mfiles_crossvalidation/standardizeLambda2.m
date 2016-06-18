function [] = standardizeLambda2()
% this function will attempt to standardize results to see if error is
% reduced with lambda2 is "low". The problem is determining what is "low".

% gathering all the data
load(sprintf('../matfiles/Xvalforcediso_LOOCV__soft_long_gauss_foriso%dkm.mat',0));
zkallh = zk_madd; zhallh = zh_Xval; ckallh = ck; vkallh = vk; 

load(sprintf('../matfiles/Xvalforcediso_LOOCV__nosoft_long_gauss_foriso%dkm.mat',0));
zkalls = zk_madd; zhalls = zh_Xval; ckalls = ck; vkalls = vk;  

zkallh = cell2mat(zkallh); zkalls = cell2mat(zkalls);
zhallh = cell2mat(zhallh); zhalls = cell2mat(zhalls);
idx = ~isnan(zkallh) & ~isnan(zkalls); 
zkallh = zkallh(idx); zhallh = zhallh(idx); 
zkalls = zkalls(idx); zhalls = zhalls(idx); 
vkallh = cell2mat(vkallh); vkallh = vkallh(idx);
vkalls = cell2mat(vkalls); vkalls = vkalls(idx);

% get standard error: (x-xhat)/s
blah = (zhallh-zkallh)./sqrt(vkallh);
blah = 5;

% coefficient of variation: mu/s















end