% expKernelSmoothTest  - Tests expKernelSmooth_stg.m & expKernelSmooth_stv.m (v.2011/01/04)

%
% Generate stg data with shift in cosinus from east to west
%
nMS=400;
nME=2*12;
Lx=5000;
Ly=2000;
rand('state',1);
randn('state',1);
xMS=Lx*rand(nMS,1);
yMS=Ly*rand(nMS,1);
cMS=[xMS,yMS];
idMS=[1:nMS]';
tME=1:nME;
TME=kron(ones(nMS,1),tME);
XMS=kron(xMS,ones(1,nME));
Z=3*randn(nMS,nME)+10*cos(2*pi*TME/12+pi*XMS/Lx);
Z=Z-quantest(Z(:),0.05);
Z(Z<0)=NaN;

%
% Space/time smoothing
%
smoothingParam=[100,50,12,3];

disp('Using expKernelSmooth_stg');
[pd,zd]=valstg2stv(Z,cMS,tME);
t1=tic;[mI]=expKernelSmooth_stg(Z,cMS,idMS,tME,smoothingParam,pd);toc(t1);
[MI]=valstv2stg(pd,mI,cMS,tME);

disp('Using expKernelSmooth_stv');
t1=tic;[mI2]=expKernelSmooth_stv(pd,zd,smoothingParam,pd);toc(t1);
[MI2]=valstv2stg(pd,mI2,cMS,tME);

for iMS=1:2, 
	figure; 
	hd=plot(tME,Z(iMS,:),'o-'); 
	hold on; 
	h1=plot(tME,MI(iMS,:),'-r'); 
	h2=plot(tME,MI(iMS,:),'--g'); 
	title(['Station ' num2str(iMS)]);
	ylabel('Air pollution concentration')
	xlabel('Time (month)')
	legend([hd h1 h2],'Data ','expKernelSmooth stg','expKernelSmooth stv');
end

