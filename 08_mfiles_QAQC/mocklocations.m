function [] = mocklocations(numlocs)
% this function will create numlocs randomly generated locations within the
% continental US

if nargin < 1, numlocs = 500; end

% axis boundaries (note these boudaries are in projected LCC format)
rand('seed',0);
ax = [ -2000000 2500000 -2000000 1500000 ];
partid = 1:numlocs; % participant ID
partx = rand(1,numlocs)*(ax(2)-ax(1)) + ax(1); % x coordinates
party = rand(1,numlocs)*(ax(4)-ax(3)) + ax(3); % y coordinates

% saves info in text file
A = [ partid ; partx ; party ];
fileID = fopen('partlocs.txt','w'); % participant locations
fprintf(fileID,'%s\t%s\t%s\n','id','x','y');
fprintf(fileID,'%d\t%0.2f\t%0.2f\n',A);
fclose(fileID);

% save in .mat file
save('../matfiles_QAQC/partlocs','partid','partx','party');

end