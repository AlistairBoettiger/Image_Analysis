%%                  view3Ddots.m
% Alistair Boettiger                            Date Begun: 02/14/11
%                                               Last Modified: 02/15/11

clear all;

folder = '/Volumes/Data/Lab Data/Raw_Data/02-06-11/';
fname = 'MP10_22C_sna_y_c';
filename = [folder,'/',fname];

Im = lsm_read_mod([filename,'.mat'],1,1.5E4); 

%%

% filtering parameters;  
alphaE = 1; %
    sigmaE = 2; %
  alphaI = 1.05; % 
   min_int  = .05; % 
   FiltSize = 25;% 
 min_size = 10; % 
    sigmaI = 3/2*sigmaE;
    Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
    Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
  

[h,w] = size(Im{1,1}{1});

z = 20; m = .85;

Il0 = Im{1,z-1}{1}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );
Il1 = Im{1,z}{1}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );
Il2 = Im{1,z+1}{1}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );

figure(1); clf; imshow(Il0);


Zs = length(Im);

[hs,ws] = size(Il1); 
%Isect = uint16(zeros(hs,ws,Zs));

Isect = zeros(hs,ws,Zs);
for z = 1:Zs % z = 20
    I =  Im{1,z}{1}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );
    
    
      % Faster method to apply filter -- use straight Gaussians.  

          outE = imfilter(single(I),Ex,'replicate'); 
          outI = imfilter(single(I),Ix,'replicate'); 
          outims = alphaE.*outE - alphaI*outI;

          bw2 = im2bw(outims,graythresh(outims)); %Automatic threshold calculated using Otsu's method.  
        %  figure(2); clf; imshow(bw2);
          D = -bwdist(~bw2);
          L = watershed(D);
          BW = bw2; BW(L==0)=0; 
         % figure(2); clf; imshow(BW);

          bw3 = im2bw(I,min_int);   
          bw2 = bw3 & BW; % Must be above threshold and shape selected by LALI filter and watershedding
         % figure(2); clf; imshow(bw3);
          
          bw2 = bwareaopen(bw2,min_size);% remove objects less than n pixels in size 
          
          bnd = imdilate(bw2,strel('disk',2)) -bw2;
          
          
         % figure(2); clf; imshow(bnd);
          
          
          mask = double(2*bw2)+bnd;   
          mask(mask==0)=NaN; 
          mask(mask==1) = 0; 
          
          % figure(2); clf; imagesc(mask);
          
    Isect(:,:,z) = double(I).*mask;
    % figure(2); clf; imagesc(  Isect(:,:,z));
    
end

[X,Y] = meshgrid((1:ws)*50,(1:hs)*50);


Istack = zeros(hs,ws,Zs);

first = 1; last = Zs;


figure(1); clf;
for z=first:last
%     I1 = double(Isect(:,:,z)); I1(I1<7000) = NaN;
%     Istack(:,:,z) = I1; 
    
    I1 = Isect(:,:,z); 
        Z = (Zs - z*ones(hs,ws))*340;
    surf(X,Y,Z,I1); hold on;
end
shading interp;
zlim([0,Zs*340]);

cmax = max(Isect(:));
C = hot(cmax);
colormap(C); colorbar; caxis([0,cmax]); set(gcf,'color','k');
view(-170,20); 
set(gca,'FontSize',12);
xlabel('nanometers');  ylabel('nanometers'); zlabel('nanometers');

axis off;

