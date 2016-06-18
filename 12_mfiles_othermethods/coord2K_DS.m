function [K]=coord2K_DS(c1,c2,model,param,param_DS,filtmodel);

% coord2K                   - covariance/variogram matrix from coordinates (Jan 1,2001)
%
% Compute the covariance or variogram matrix between two sets of coordinates,
% based on the Euclidean distances between these sets of coordinates.
%
% SYNTAX :
%
% [K]=coord2K(c1,c2,model,param,filtmodel);
%
% INPUT :
%
% c1        n1 by d  matrix of coordinates for the locations in the first set.
%                    A line corresponds to the vector of coordinates at a location,
%                    so the number of columns is equal to the dimension of the space.
%                    There is no restriction on the dimension of the space.
% c2        n2 by d  matrix of coordinates for the locations in the second set, using
%                    the same conventions as for c1. 
% model     string   that contains the name of the variogram or covariance model that
%                    is used for the computation.
% param     k by 1   vector of values for the parameters of model, according to the
%                    convention for the corresponding variogram or covariance model
%                    (see the MODELS directory).
% filtmodel scalar   optional value specifying if the stochastic model component is to
%                    be included (filtmodel=1) or is to be filtered out (filtmodel=0).
%                    See the krigingfilter.m function in the KRIGING directory.
%                    
%
% OUTPUT :
%
% K         n1 by n2 covariance or variogram matrix between the coordinates in c1 and c2,
%                    depending on the fact that model is a covariance or a variogram model.
%                    The number of lines and columns for K corresponds to the number of
%                    lines in c1 and the number of lines in c2, respectively.
%
% NOTE :
%
% For a detailed discussion about the coding of the model and param variables in
% various situations (e.g., nested models, multivariate case, space/time case),
% the reader is referred to the detailed description given for the kriging.m
% function in the KRIGING directory.

%%%%%% Build the general data structure for singular cases

noindex=(iscell(c1)==0);
nocell=(iscell(model)==0);

if nocell==1 & noindex==1,          % case for 1 model and 1 variable
  n1=size(c1,1);                    % compute the number of coordinates
  n2=size(c2,1);
  c1={c1,ones(n1,1)};               % create an index with values 1
  c2={c2,ones(n2,1)};
  nm=1;                             % set number of models equal to 1
  nv=1;                             % set number of variables equal to 1
  model={model};                    % create a cell array with model
  np=length(param);                 % create a cell array with param
  param={{param(1),param(2:np)}};
end;

if nocell==1 & noindex==0,          % case for 1 model and nv variables
  n1=size(c1{1},1);                 % compute the number of coordinates
  n2=size(c2{1},1);
  nm=1;                             % set number of models equal to 1
  nv=max([c1{2};c2{2}]);            % compute the number of variables 
  model={model};                    % create a cell array with model
  param={param};
end;

if nocell==0 & noindex==1,          % case for nm models and 1 variable
  n1=size(c1,1);                    % compute the number of coordinates
  n2=size(c2,1);
  c1={c1,ones(n1,1)};               % create an index with values 1
  c2={c2,ones(n2,1)};
  nm=length(model);                 % compute the number of models 
  nv=1;                             % set the number of variables to 1 
  for i=1:nm,                       % create cell arrays with param vectors 
    np=length(param{i});               
    param{i}={param{i}(1),param{i}(2:np)};
  end;
end;

if nocell==0 & noindex==0;          % case for nm models and nv variables
  n1=size(c1{1},1);                 % compute the number of coordinates
  n2=size(c2{1},1);
  nm=length(model);                 % compute the number of models
  nv=max([c1{2};c2{2}]);            % compute the number of variables 
end;

%%%%%% Initialize the parameters

if nargin<6,
  filtmodel=ones(nm,1);
end;

%%%%%% compute the distance matrix and initialize the size of K

[isST,isSTsep,modelS,modelT]=isspacetime(model);

if isSTsep,
  [nparam,models]=nparammodels;
end;

if (n1==0)|(n2==0),
  K=zeros(n1,n2);
  return;
else
  if ~isST,
    D=coord2dist(c1{1},c2{1});
    K=zeros(size(D));
  else
    nd=size(c1{1},2)-1;
    Ds=coord2dist(c1{1}(:,1:nd),c2{1}(:,1:nd));
    Dt=coord2dist(c1{1}(:,nd+1),c2{1}(:,nd+1));
    K=zeros(size(Ds));
  end;
end;

%%%%%% Compute the covariance matrix or vector

if n1==n2,                            % test if K is symmetric and set
  issymm=prod(prod(double(c1{1}==c2{1})));    % issymm to 0 or 1 
else
  issymm=0;
end;

if issymm==1,                 % if K is symmetric 
  for i=1:nv,                 % loop twice over the number of variables
    indexi=find(c1{2}==i);    % and the number of models by selecting the 
    for j=i:nv,               % appropriate subset of locations and model
      indexj=find(c2{2}==j);  % and use the symmetry property of the matrix
      for k=1:nm,
        if filtmodel(k)==1,
	  C=param_DS{k};
          if ~isST,
              if n1 == 1 & n2 == 1 & ~isempty(param{k}{2}) % variance at estimation location
                  K(indexi,indexj)=K(indexi,indexj)+eval([model{k},'(D(indexi,indexj),{C(indexi,indexj),param{k}{2}})']);
              elseif length(C) == 1 % old way
                K(indexi,indexj)=K(indexi,indexj)+eval([model{k},'(D(indexi,indexj),[C,param{k}{2}])']);
              else % the Cd by Cd calculation
                K(indexi,indexj)=K(indexi,indexj)+eval([model{k},'(D(indexi,indexj),{C(indexi,indexj),param{k}{2}})']);
              end
          end;
          if (~isSTsep)&(isST),
	    K(indexi,indexj)=K(indexi,indexj)+eval...
            ([model{k},'(Ds(indexi,indexj),Dt(indexi,indexj),[C(i,j),param{k}{2}])']);
          end;
          if isSTsep,              
              if length(C) == 1 % old way
                nps=nparam{strmatch(modelS{k},models,'exact')};
                npt=nparam{strmatch(modelT{k},models,'exact')};	    
                Ks=eval([modelS{k},'(Ds(indexi,indexj),[1,param{k}{2}(1:nps-1)])']); % old way
                Kt=eval([modelT{k},'(Dt(indexi,indexj),[1,param{k}{2}(nps:nps+npt-2)])']); % old way
                K(indexi,indexj)=K(indexi,indexj)+C(i,j)*Ks.*Kt; % old way
              elseif n1 == 1 & n2 == 1 & ~isempty(param{k}{2}) % variance at estimation location
                Ks=eval([modelS{k},'(Ds(indexi,indexj),{1,param{k}{2}(1)})']); % new way
                Kt=eval([modelT{k},'(Dt(indexi,indexj),[1,param{k}{2}(2)])']); % new way
                K(indexi,indexj)=K(indexi,indexj)+C(indexi,indexj).*Ks.*Kt; % new way
              else % the Cd by Cd calculation
                Ks=eval([modelS{k},'(Ds(indexi,indexj),{1,param{k}{2}(1)})']); % new way
                Kt=eval([modelT{k},'(Dt(indexi,indexj),[1,param{k}{2}(2)])']); % new way
                K(indexi,indexj)=K(indexi,indexj)+C(indexi,indexj).*Ks.*Kt; % new way
              end
          end;
        end;
      end;
      if i~=j,
        K(indexj,indexi)=K(indexi,indexj)';
      end;
    end;
  end;
else                          % else K is not symmetric
  for i=1:nv,                 % loop twice over the number of variables
    indexi=find(c1{2}==i);    % and the number of models by selecting the
    for j=1:nv,               % appropriate subset of locations and model
      indexj=find(c2{2}==j);
      for k=1:nm,
        if filtmodel(k)==1,
	  C=param_DS{k};
          if ~isST,
              if length(C) == 1 % old way
                K(indexi,indexj)=K(indexi,indexj)+eval([model{k},'(D(indexi,indexj),[C,param{k}{2}])']);
              elseif size(C,1) ~= size(C,2) % Ck by Cd calculation 
                  K(indexi,indexj)=K(indexi,indexj)+eval([model{k},'(D(indexi,indexj),{C(indexj,indexi)'',param{k}{2}})']);
              else % Cd by Cd calculation
                a = combvec(indexi',indexj'); x = a(2,:)'; y = a(1,:)'; 
                temp1 = cell2mat( arrayfun(@(x,y) eval([model{k},'(D(x,y),[C(x,y),param{k}{2}])']),x,y,'UniformOutput',false) );
                temp2 = reshape(temp1,length(indexi),length(indexj));
                K(indexi,indexj)=K(indexi,indexj)+temp2;
              end            
          end;
          if (~isSTsep)&(isST),
	    K(indexi,indexj)=K(indexi,indexj)+eval...
            ([model{k},'(Ds(indexi,indexj),Dt(indexi,indexj),[C(i,j),param{k}{2}])']);
          end;
          if isSTsep,
              if length(C) == 1 % old way
                nps=nparam{strmatch(modelS{k},models,'exact')};
                npt=nparam{strmatch(modelT{k},models,'exact')};
                Ks=eval([modelS{k},'(Ds(indexi,indexj),[1,param{k}{2}(1:nps-1)])']);
                Kt=eval([modelT{k},'(Dt(indexi,indexj),[1,param{k}{2}(nps:nps+npt-2)])']);
                K(indexi,indexj)=K(indexi,indexj)+C(i,j)*Ks.*Kt;
              elseif size(C,1) ~= size(C,2) % Ck by Cd calculation 
                Ks=eval([modelS{k},'(Ds(indexi,indexj),{1,param{k}{2}(1)})']); % new way
                Kt=eval([modelT{k},'(Dt(indexi,indexj),[1,param{k}{2}(2)])']); % new way
                K(indexi,indexj)=K(indexi,indexj)+C(indexj,indexi)'.*Ks.*Kt; % new way             
              else % the Cd by Cd calculation
                a = combvec(indexi',indexj'); x = a(2,:)'; y = a(1,:)';
                Ks=eval([modelS{k},'(Ds(indexi,indexj),{1,param{k}{2}(1)})']); % new way
                Kt=eval([modelT{k},'(Dt(indexi,indexj),[1,param{k}{2}(2)])']); % new way
                temp1 = cell2mat( arrayfun(@(x,y) C(x,y),x,y,'UniformOutput',false) );
                temp2 = reshape(temp1,length(indexi),length(indexj));
                K(indexi,indexj)=K(indexi,indexj)+temp2.*Ks.*Kt; % new way
              end
          end;
        end;
      end;
    end;
  end;
end;

