function plotmask(pVal,ax,maskfillcolor,masklinecolor)
% plotMask      - Mask the map
%
% SYNTAX :
% 
%  plotmask(pVal,ax,maskfillcolor,masklinecolor)
% 
% INPUT :
% 
%  pVal                 :Cordinate values read from e00 file using readARCe00 function
%  ax                   :The area you interested in
%                        Specify the four elements vector [minLongitude maxLongitude minLatitude maxLatitude]
%  maskfillcolor        :Maskfillcolor, Default is white
%  masklinecolor        :Masklinecolor, Default is black
%
% OUTPUT :
%
%  maskdata.mat         :The file is created in the maskfile folder
%                        This file contains following parameters
%                        You can change the name and the path of this file
%                        These are set in Line 40
%
%   'maskedge'          :This parameter used to draw the most outer edge of the mask
%   'idxIns'            :This parameter used to draw the area inside the main area
%   'pArea'             :This parameter used to draw the outer edge of each area
%   'ax'                :This parameter is the same as specified in INPUT section
%   'nArea'             :This parameter used to draw the outer edge of each area
%
%  sample code to redraw the mask using the output file is discribed below
%
%
%%%%(S A M P L E  C O D E)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% load(maskdata.mat)
%
% maskfillcolor='w';
% masklinecolor='k';
%
% figure;
% hold on;
%
% fill(maskedge(:,1),maskedge(:,2),maskfillcolor,'EdgeColor',maskfillcolor);
% 
% if length(idxIns)~=0
%     for i=1:length(idxIns)
%         fill(pArea{idxIns(i)}(:,1),pArea{idxIns(i)}(:,2),maskfillcolor,'EdgeColor',maskfillcolor);
%     end
% end
% 
% for i=1:nArea
%     plot(pArea{i}(:,1),pArea{i}(:,2),masklinecolor);
% end
% 
% axis equal;
% axis(ax);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if nargin<1, error('You have to specify the pVal'); end;
if nargin<2, error('You have to specify ax'); end;
if nargin<3, maskfillcolor='w'; end;
if nargin<4, masklinecolor='k'; end;

maskedgefile=['maskdata'];

nArea=1;
pArea{nArea}=[pVal{1}];

for i=1:length(pVal)-1
    if pVal{i}(end,1)==pVal{i+1}(1,1) & pVal{i}(end,2)==pVal{i+1}(1,2)
        pArea{nArea}=[pArea{nArea};pVal{i+1}];
    else
        nArea=nArea+1;
        pArea{nArea}=[pVal{i+1}];
    end
end

for i=1:nArea
    idxsameval=[find(sum(abs(diff(pArea{i})),2)~=0);size(pArea{i},1)];
    pArea{i}=pArea{i}(idxsameval,:);
    minAreaLo(i)=min(pArea{i}(:,1));
    maxAreaLo(i)=max(pArea{i}(:,1));
    minAreaLa(i)=min(pArea{i}(:,2));
    maxAreaLa(i)=max(pArea{i}(:,2));
    idxstpt(i)=min(find(pArea{i}(:,1)==min(pArea{i}(:,1))));
end

puax=ax;

if ax(1)>min(minAreaLo)
    puax(1)=min(minAreaLo)-(max(maxAreaLo)-min(minAreaLo))/20;
end

if ax(2)<max(maxAreaLo)
    puax(2)=max(maxAreaLo)+(max(maxAreaLo)-min(minAreaLo))/20;
end

if ax(3)>min(minAreaLa)
    puax(3)=min(minAreaLa)-(max(maxAreaLo)-min(minAreaLo))/20;
end

if ax(4)<max(maxAreaLa)
    puax(4)=max(maxAreaLa)+(max(maxAreaLo)-min(minAreaLo))/20;
end

[srtLo,idxLo]=sort(minAreaLo);
maskedge=[puax(1) puax(4);puax(1) pArea{idxLo(1)}(idxstpt(idxLo(1)),2)];

idxinArea=[];
for i=1:nArea
    for j=1:nArea
        if i==j
            idxinArea(i,j)=0;
        else
            if maxAreaLo(idxLo(i))>maxAreaLo(idxLo(j))
                flgin=sum(inpolygon(pArea{idxLo(j)}(:,1),pArea{idxLo(j)}(:,2),...
                    pArea{idxLo(i)}(:,1),pArea{idxLo(i)}(:,2)));
                if flgin==0
                    idxinArea(i,j)=0;
                else
                    idxinArea(i,j)=1;
                end
            else
                idxinArea(i,j)=0;
            end
        end
    end
end
idxinside=sum(idxinArea);

idxSep=idxLo(~logical(idxinside));
idxIns=idxLo(logical(idxinside));

if length(idxSep)==1
    if idxstpt(idxSep(1))~=1
        maskedge=[maskedge;pArea{idxSep(1)}(idxstpt(idxSep(1)):end,:);...
                pArea{idxSep(1)}(2:idxstpt(idxSep(1)),:)];
    else
        maskedge=[maskedge;pArea{idxSep(1)}(:,:)];
    end
    maskedge=[maskedge;puax(1) pArea{idxSep(1)}(idxstpt(idxSep(1)),2);...
            puax(1) puax(3);puax(2) puax(3);puax(2) puax(4);puax(1) puax(4)];
elseif length(idxSep)==2
    len2=[];
    for i=1:length(pArea{idxSep(1)})
        len1=[];
        for j=1:length(pArea{idxSep(2)})
            len1(j)=sqrt((pArea{idxSep(1)}(i,1)-pArea{idxSep(2)}(j,1))^2+...
                    (pArea{idxSep(1)}(i,2)-pArea{idxSep(2)}(j,2))^2);
        end
        [minlen1,idxlen1]=min(len1);
        len2(i,:)=[minlen1 idxlen1];
    end
    [minlen2,idxlen2]=min(len2(:,1));
    pton1=idxlen2;
    pton2=len2(idxlen2,2);
    
    if idxstpt(idxSep(1))<pton1
        maskedge=[maskedge;pArea{idxSep(1)}(idxstpt(idxSep(1)):pton1,:)];
    else
        maskedge=[maskedge;pArea{idxSep(1)}(idxstpt(idxSep(1)):end,:);...
                pArea{idxSep(1)}(2:pton1,:)];
    end
    
    if pton2==1
        maskedge=[maskedge;pArea{idxSep(2)}(pton2:end,:)];
    else
        maskedge=[maskedge;pArea{idxSep(2)}(pton2:end,:);pArea{idxSep(2)}(2:pton2,:)];
    end
    
    if idxstpt(idxSep(1))<pton1
        maskedge=[maskedge;pArea{idxSep(1)}(pton1:end,:);...
                pArea{idxSep(1)}(2:idxstpt(idxSep(1)),:)];
    else
        maskedge=[maskedge;pArea{idxSep(1)}(pton1:idxstpt(idxSep(1)),:)];
    end

    maskedge=[maskedge;puax(1) pArea{idxSep(1)}(idxstpt(idxSep(1)),2);...
            puax(1) puax(3);puax(2) puax(3);puax(2) puax(4);puax(1) puax(4)];
else
    for i=1:length(idxSep)-1
        if maxAreaLo(idxSep(i))>minAreaLo(idxSep(i+1))
            error('This function cannot deal with such a complicated shape')
        end
        len2=[];
        for j=1:length(pArea{idxSep(i)})
            len1=[];
            for k=1:length(pArea{idxSep(i+1)})
                len1(k)=sqrt((pArea{idxSep(i)}(j,1)-pArea{idxSep(i+1)}(k,1))^2+...
                    (pArea{idxSep(i)}(j,2)-pArea{idxSep(i+1)}(k,2))^2);
            end
            [minlen1,idxlen1]=min(len1);
            len2(j,:)=[minlen1 idxlen1];
        end
        [minlen2,idxlen2]=min(len2(:,1));
        len3(i,:)=[minlen2,idxlen2,len2(idxlen2,2)];
    end

    if idxstpt(idxSep(1))<len3(1,2)
        maskedge=[maskedge;pArea{idxSep(1)}(idxstpt(idxSep(1)):len3(1,2),:)];
    else
        maskedge=[maskedge;pArea{idxSep(1)}(idxstpt(idxSep(1)):end,:);...
                pArea{idxSep(1)}(2:len3(1,2),:)];
    end
    
    for i=2:length(idxSep)-1
        if len3(i-1,3)<len3(i,2)
            maskedge=[maskedge;pArea{idxSep(i)}(len3(i-1,3):len3(i,2),:)];
        else
            maskedge=[maskedge;pArea{idxSep(i)}(len3(i-1,3):end,:);...
                    pArea{idxSep(i)}(2:len3(i,2),:)];
        end
    end
    
    if len3(end,3)==1
        maskedge=[maskedge;pArea{idxSep(end)}(len3(end,3):end,:)];
    else
        maskedge=[maskedge;pArea{idxSep(end)}(len3(end,3):end,:);pArea{idxSep(end)}(2:len3(end,3),:)];
    end
    
    for i=length(idxSep)-1:-1:2
        if len3(i-1,3)<len3(i,2)
            maskedge=[maskedge;pArea{idxSep(i)}(len3(i,2):end,:);...
                    pArea{idxSep(i)}(2:len3(i-1,3),:)];
        else
            maskedge=[maskedge;pArea{idxSep(i)}(len3(i,2):len3(i-1,3),:)];
        end
    end
    
    if idxstpt(idxSep(1))<len3(1,2)
        maskedge=[maskedge;pArea{idxSep(1)}(len3(1,2):end,:);...
                pArea{idxSep(1)}(2:idxstpt(idxSep(1)),:)];
    else
        maskedge=[maskedge;pArea{idxSep(1)}(len3(1,2):idxstpt(idxSep(1)),:)];
    end

    maskedge=[maskedge;puax(1) pArea{idxSep(1)}(idxstpt(idxSep(1)),2);...
            puax(1) puax(3);puax(2) puax(3);puax(2) puax(4);puax(1) puax(4)];
end

% save(maskedgefile,'maskedge','idxIns','pArea','ax','nArea');

hold on;
fill(maskedge(:,1),maskedge(:,2),maskfillcolor,'EdgeColor',maskfillcolor);

if length(idxIns)~=0
    for i=1:length(idxIns)
        fill(pArea{idxIns(i)}(:,1),pArea{idxIns(i)}(:,2),maskfillcolor,'EdgeColor',maskfillcolor);
    end
end

for i=1:nArea
    plot(pArea{i}(:,1),pArea{i}(:,2),masklinecolor);
end

%axis equal;
axis(ax);