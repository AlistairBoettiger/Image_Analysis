%%                  view3Ddots.m
% Alistair Boettiger                            Date Begun: 02/14/11
%                                               Last Modified: 02/15/11

clear all;

% folder = '/Volumes/Data/Lab Data/Raw_Data/02-06-11/';
% fname = 'MP10_22C_sna_y_c'; emb = '01';

folder = '/Volumes/Data/Lab Data/Raw_Data/02-17-11/';
fname = 'MP05_22C_sna_y'; emb = '02';
%fname = 'MP09_22C_hb_y_d'; emb = '01';

filename = [folder,'/',fname];

Im = lsm_read_mod([filename,'.mat'],str2double(emb),1.5E4); 

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

z = 20; m = .88; % .85;

Il0 = Im{1,z-1}{1}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );
Il1 = Im{1,z}{1}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );
Il2 = Im{1,z+1}{1}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );

figure(1); clf; imshow(Il0);


Zs = length(Im);

[hs,ws] = size(Il1); 
%Isect = uint16(zeros(hs,ws,Zs));

Isect1 = zeros(hs,ws,Zs);
Isect2 = zeros(hs,ws,Zs);
Inuc = zeros(hs,ws,Zs);

for z = 1:Zs % z = 20        
      % Dot finding for channel 1
          I1 =  Im{1,z}{1}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );
          [cent1,bw1] = dotfinder(I1,alphaE,alphaI,Ex,Ix,min_int,min_size);
          bnd1 = imdilate(bw1,strel('disk',2)) -bw1;         
         % figure(2); clf; imshow(bnd2);
          mask = double(2*bw1)+bnd1;   
          mask(mask==0)=NaN; 
          mask(mask==1) = 0;       
          % figure(2); clf; imagesc(mask);          
          Isect1(:,:,z) = 255*double(I1)/2^16.*mask;
          % figure(2); clf; imagesc(  Isect1(:,:,z));
         
          % Dot finding for channel 2
          I2 = Im{1,z}{2}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) );
          [cent2,bw2] = dotfinder(I2,alphaE,alphaI,Ex,Ix,min_int,min_size);
          bnd2 = imdilate(bw2,strel('disk',2)) -bw2;         
         % figure(2); clf; imshow(bnd2);
          mask = double(2*bw2)+bnd2;   
          mask(mask==0)=NaN; 
          mask(mask==1) = 0;       
          % figure(2); clf; imagesc(mask);          
          Isect2(:,:,z) = 255*double(I2)/2^16.*mask;
          
    % Nuclear finding            
    In = Im{1,z}{3}( floor(h/2*m):floor(h/2*(2-m)), floor(w/2*m):floor(w/2*(2-m)) ); 
    nucnoise = double(im2bw(In,.15)); 
    nucnoise(nucnoise==0) = NaN; 
    Inuc(:,:,z) = 255*double(In)/2^16.*nucnoise;
    figure(2); clf; imagesc(Inuc(:,:,z));
    
    
end


%%
[X,Y] = meshgrid((1:ws)*50,(1:hs)*50);
Istack = zeros(hs,ws,Zs);
first = 1; last = Zs;

nmax = round(max(Inuc(:)));
c1max = round(max(Isect1(:)));
c2max = round(max(Isect2(:)));

% 
figure(1);  clf;
colordef black; set(gca,'color','k'); set(gcf,'color','k');

% Plot nuclei data
for z=first:last
    In = Inuc(:,:,z) + c1max + c2max; 
    Z = (Zs - z*ones(hs,ws))*340;
    surf(X,Y,Z,In); hold on;
end
shading interp;
alpha(.45);  % make nuclei transparent 



figure(1);  hold on;
for z=first:last
    I1 = Isect1(:,:,z) ;
   Z = (Zs - z*ones(hs,ws))*340;
    surf(X,Y,Z,I1); hold on;
end
shading interp;

figure(1);  hold on;
for z=first:last
    I2 = Isect2(:,:,z) + c1max+1; 
    Z = (Zs - z*ones(hs,ws))*340;
    surf(X,Y,Z,I2); hold on;
end
shading interp;

C2 = zeros(c2max,3);
for cc = 1:c2max
    C2(cc,:) = [((cc-1)/c2max)^3,((cc-1)/c2max)^.5,((cc-1)/c2max)^2];
end

CN = cool(nmax); 
C1 = hot(c1max);
%C1 = flipud(1-hot(c1max));


colormap([0,0,0;C1;C2;CN]); colorbar; caxis([0,nmax+c1max+c2max]); 
%view(-109,50); 
view(40,39); axis on;
set(gca,'FontSize',12);
zlim([0,Zs*340]);
xlabel('nanometers');  ylabel('nanometers'); zlabel('nanometers');

axis off;

