% smoothingExample  - Example of using expKernelSmooth_stg vs stmean_pI


stmean_method='stmean_pI'; % Use stmean_method='stmean' to use the old 
                           %   stmean function, which is quicker
												   % Use stmean_method='stmean_pI' to use the new  
												   %   stmean function, which is a little more
													 %   accurate but slower

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
smoothingParam=[300,50,12,3];

disp('Using expKernelSmooth_stg');
[pd,zd]=valstg2stv(Z,cMS,tME);
t1=tic;[mI1]=expKernelSmooth_stg(Z,cMS,idMS,tME,smoothingParam,pd);toc(t1);
[MI1]=valstv2stg(pd,mI1,cMS,tME);

disp(stmean_method);
switch stmean_method
	case 'stmean'
		[ms,mss,mt,mts]=stmean(Z,cMS,idMS,tME,smoothingParam);
		t1=tic;[mI2]=stmeaninterpstv(cMS,tME,mss,mts,pd);toc(t1);
	case 'stmean_pI'
		t1=tic;[mI2]=stmean_pI(Z,cMS,idMS,tME,smoothingParam,pd);toc(t1);
	otherwise
		error('stmean_method must be set to ''stmean'' or ''stmean_pI''');
end
[MI2]=valstv2stg(pd,mI2,cMS,tME);

for iMS=1:10, 
	figure;
	hd=plot(tME,Z(iMS,:),'o-');
	hold on;
	hI1=plot(tME,MI1(iMS,:),'-r');
	hI2=plot(tME,MI2(iMS,:),':k');
	title(['Station ' num2str(iMS)]);
	ylabel('Air pollution concentration')
	xlabel('Time (month)')
	legend([hd hI1 hI2],['Data at station ' num2str(iMS)],'kernel smoothing','stmean');
	print('-dpng',sprintf('station%d.png',iMS));
end






