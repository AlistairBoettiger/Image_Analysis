
function plotstackdots(Isect1,varargin)
%-------------------------------------------------------------------------%
% Output
%        figure(1) shows 3D reconstruction of confocal sections;
%-------------------------------------------------------------------------%
% Input 
%    Isect1 -- HxWxZs matrix, where H and W are the dimensions of the image
%               and Zs is the number of z-sections in the confocal stack.
%               It is best that H and W are small (i.e. subset of image).
%               zeros in this matrix will be invisible.  nonzeros will
%               appear as solid sections.  Set figure(1) alpha < 1 to make
%               nonzero values semi-transparent. 
%               This variable is inteded for mRNA spot data. 
%-------------------------------------------------------------------------%
% Optional Inputs
%    Isect2 -- HxWxZs matrix, just like Isect1.  Will be rendered green in
%               the reconstructed image
%    Inuc --   HxWxZs matrix, just like Isect1.  Alpha will be set to .45
%               by default and color to blue.  (This matrix is generally
%               used for nuclei data).  
%-------------------------------------------------------------------------%
%
% Alistair Boettiger
% Version 1.0
% Copyright Creative Commons 3.0 CC BY  September 26, 2012.
% 
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Parse variable inupts
%-------------------------------------------------------------------------%
if nargin == 1
    Isect2 = NaN*ones( size(Isect1) );
    Inuc =  NaN*ones(size(Isect1) );
elseif nargin == 2
    Isect2 = varargin{1};
    Inuc =  NaN*ones(size(Isect1) );
elseif nargin ==3
    Isect2 = varargin{1};
    Inuc = varargin{2};
else
    error('expected 1-3 doubles of size H x W x Zs: mRNA 1, mRNA 2, and Nnucs'); 
end

if isempty(Isect2)
    Isect2 =  NaN*ones( size(Isect1) );
end
%-------------------------------------------------------------------------%
%% Main Function
%-------------------------------------------------------------------------%

% Ensure the correct class (otherwise NaNs become 0s)
Isect1 = double(Isect1);
Isect2 = double(Isect2);
Inuc = double(Inuc);

% leave blank regions not containing dots
Isect1(Isect1==0) = NaN;
Isect2(Isect2==0) = NaN;
Inuc(Inuc==0) = NaN; 
[hs,ws,Zs] = size(Isect1); 


% Build grid
[X,Y] = meshgrid((1:ws)*50,(1:hs)*50);
first = 1; last = Zs;

% to help specify colors correctly
nmax = max(0,round(max(Inuc(:))));
c1max = max(0,round(max(Isect1(:))));
c2max = max(0,round(max(Isect2(:))));
% max(0,x) gets around NaNs

% 
figure(1);  clf;
colordef black; 
set(gca,'color','k');
set(gcf,'color','k');

% Plot nuclei data
for z=first:last
    In = Inuc(:,:,z) + c1max + c2max; 
    Z = (Zs - z*ones(hs,ws))*340;
    surf(X,Y,Z,In); hold on;
end
shading interp;
 alpha(.45);  % make nuclei transparent 


% Plot red channel confocal data
figure(1);  hold on;
for z=first:last
    I1 = Isect1(:,:,z) ;
   Z = (Zs - z*ones(hs,ws))*340;
    surf(X,Y,Z,I1); hold on;
end
shading interp;

% Plot green channel confocal data
figure(1);  hold on;
for z=first:last
    I2 = Isect2(:,:,z) + c1max+1; 
    Z = (Zs - z*ones(hs,ws))*380;
    surf(X,Y,Z,I2); hold on;
end
shading interp;

% build colormaps
C2 = zeros(c2max,3);
for cc = 1:c2max
    C2(cc,:) = [((cc-1)/c2max)^3,((cc-1)/c2max)^.5,((cc-1)/c2max)^2];
end
CN = cool(nmax); 
C1 = hot(c1max);
colormap([0,0,0;C1;C2;CN]); 
colorbar; 
caxis([0,nmax+c1max+c2max]); 

view(-109,50); 
 axis off;
