
% Alistair Boettiger                                Date Begun: 09/25/12
% Levine Lab            

%% Attempt for cover: Cell Reports paper
%
%% Publication figures as of 12/14/12
clear all;


% % im 2 (larger region for nuclei stacks)
imdate = '2011-02-17/'; stackfolder = 'MP05_22C\';
fname = 'MP05_22C_sna_y_c'; rawfolder = 'G:\Raw_Data\';
emb = '02'; 
    sigmaE =  2.5; % 3;%   IMPORTANT    3 for LSM700, 2.5 for LSM710
    sigmaI =  3.5; % 4; %  IMPORTANT
    min_int  = 0.015;    %  5    ;% .05 % not necessary Fix at Zero
    FiltSize = 30;% 
    min_size = 20; % 30;% 
    min_peak = 1000;%
    max_size = 110; % 130; 
folder = 'C:\Users\Alistair\My Documents\Projects\mRNA_counting\Data\'; 
rawfolder = [rawfolder,'/',imdate];
folder = [folder,'/',imdate];
load([folder,'/',fname,'_',emb,'_nucdata.mat']);
h1 = 1200;
h2 = 1600;
w1 = 200;
w2 = 600;
hs = 400;% 150;
ws = 400;% 150;
hz1 = 50;
hz2 = hz1+hs;
wz1 = 50;
wz2 = wz1+ws;
%     % ------

Zs = 35;
addpath('C:\Users\Alistair\Documents\Projects\mRNA_counting\Code\');
addpath('C:\Users\Alistair\Documents\Projects\Image_Analysis\');
Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
Isect1 = zeros(hs,ws,Zs); 
Inuc = zeros(hs,ws,Zs); 
im_folder = cell(Zs,1); 

for z= 1:Zs;  % 14 10
    %%
im_folder{z} = [rawfolder,stackfolder,fname,'_',emb,'_z',num2str(z),'.tif'];
Iin_z = imreadfast(im_folder{z});               


I3 = Iin_z(h1:h2,w1:w2,:);


% A: Original image
I3(:,:,1) = mycontrast(I3(:,:,1),.0001,.1);
I3(:,:,2) = 0*I3(:,:,2);
 I3(:,:,3) = mycontrast(I3(:,:,3),.001,.01);
figure(1); clf; colordef black;  set(gcf,'color','k'); imagesc(I3); 

% array
figure(10); clf; subplot(3,2,1); 
colordef black;  set(gcf,'color','k'); imagesc(I3); 
title(['single raw image section z=',num2str(z)]);

%......
  I =  Iin_z((h1+1):h2,(w1+1):w2,1);
  bw = im2bw(I,min_int);  % if we're gonna do this, do it first.   
  Iin = I;
  Iin(bw==0) =0; 
  clear bw;  
  outE = imfilter(single(Iin),Ex,'replicate'); 
  outI = imfilter(single(Iin),Ix,'replicate'); 
  Iin = makeuint(outE - outI,16);
%.......

% B: Difference of gaussian filter
figure(2); clf; colordef black;
 It = outE - outI;
 imagesc(It); 
 shading flat;  colorbar; 
 colormap hot; caxis([.08*min(It(:)),max(It(:))]);
 set(gcf,'color','k');     
 
 figure(10); subplot(3,2,2);  imagesc(It); 
 shading flat;  colormap hot;
 caxis([.08*min(It(:)),max(It(:))]);
 set(gcf,'color','k');  
 title('Difference of Gaussians filter');  
%..........
% automatic thresholding
    [h,w] = size(I); 
    if h>400
        m=.8;
    else
        m=1/h;
    end
    xp1= floor(h/2*m)+1; 
    xp2 = floor(h/2*(2-m))+1;
    yp1 = floor(w/2*m)+1;
    yp2 = floor(w/2*(2-m))+1;

    N = 100; 
    count = zeros(1,N); 
    thresh = linspace(0,1,N);
    for t=1:N
        bw = im2bw(Iin(xp1:xp2,yp1:yp2),thresh(t)); 
        [jnk,count(t)] = bwlabel(bw); 
    end
    [jnk,im] = max(count);
    imthresh = thresh(im);
    bw2 = im2bw(Iin,imthresh);
    Iin = Iin.*uint16(bw2);  
 %..........
 
% C: Automatic threshold determination
figure(3); clf; colordef black;
plot(thresh,count,'linewidth',2); 
set(gcf,'color','k'); xlabel('theta'); 
ylabel('Number of objects');  
title('automatic threshold computation');

figure(10); subplot(3,2,3); 
plot(thresh,count,'linewidth',2); 
set(gcf,'color','k'); xlabel('theta'); 
ylabel('Number of objects');  
title('automatic threshold computation');

% D: Results of thresholding
figure(4); clf; colordef black;
imagesc(Iin); shading flat; 
colormap hot; set(gcf,'color','k'); 
title('thresholded image');

figure(10); subplot(3,2,4); 
imagesc(Iin); shading flat; 
colormap hot; set(gcf,'color','k'); 
title('thresholded image');
 
%.............   
% watershed separation
 M = max(Iin(:)); 
 L = watershed(double(M-Iin));  %  figure(2); clf; imagesc(L); colormap lines; shading flat;
 Iin(L==0) = 0;
 bw2 = logical(Iin); 
 bw2 = bwareaopen(bw2,min_size);% remove objects less than n pixels in size 
 Iin = Iin.*uint16(bw2);
%.............

% E: Watershed
figure(5); clf; colordef black;
imagesc(makeuint(Iin,16)); shading flat; 
colormap hot; set(gcf,'color','k'); 
title('watershed filter');

figure(10); subplot(3,2,5); 
imagesc(Iin); shading flat; 
colormap hot; set(gcf,'color','k'); 
title('watershed filter');
% ............. 

   bw2 = logical(Iin); 
   bw2 = bwareaopen(bw2,min_size);% remove objects less than n pixels in size 
   labeled = bwlabel(bw2,8); % count and label RNAs (8-> diagnols count as touch)  
   R1 = regionprops(labeled,Iin,'WeightedCentroid','Area'); % compute mRNA centroids      
   cent = reshape([R1.WeightedCentroid],2,length(R1))';
   
   fckups = find([R1.Area]>max_size);
   for i=1:length(fckups)
       bw2(labeled==fckups(i)) = 0;
   end
   Iin = Iin.*uint16(bw2);
   figure(4); clf; imagesc(Iin);
   
       inds = floor(cent(:,2))+floor(cent(:,1))*h;  % indices in this layer  
       inds(inds>w*h) = w*h; 
       ints = [I(inds)]; 
     
       good_dots_index = ints>min_peak;
       cent = cent(good_dots_index,:);  % EXPORTED VAR 
       inds = [inds(good_dots_index)]; % EXPORTED VAR 
       ints = [0; ints(good_dots_index)];  % EXPORTED VAR
       dot_labels = 1:length(inds); % EXPORTED VAR  
% ................     
     
% panelF = I3;
% panelF(:,:,1) = panelF(:,:,1) + uint16(2^14*cell_bounds);
% panelF(:,:,2) = panelF(:,:,2) + uint16(2^14*cell_bounds);
% panelF(:,:,3) = panelF(:,:,3) + uint16(2^14*cell_bounds);
% 
% % F: mRNA positions and nuclear assignments
%       figure(10); subplot(3,2,6);
%       imagesc(panelF); hold on;
%       plot(cent(:,1),cent(:,2),'w+');
%       title('mRNA positions');
      
 Isect1(:,:,z) = Iin; % (hz1:hz2,wz1:wz2);  %  Iin(hz1+1:hz2,wz1+1:wz2); 

 In = Iin_z(hz1+1:hz2,wz1+1:wz2,3); 
 % figure(3); clf; imagesc(In);
 nmask = im2bw(In,.2);
 nmask = imdilate(nmask,strel('disk',2));
 nmask = bwareaopen(nmask,100);
 figure(3); clf; imagesc(nmask);
 In(In<median(single(In(:)))) = 0;
 Inuc(:,:,z) = In.*uint16(nmask); 
 %%
end      

plotstackdots(Isect1,[],Inuc)
view(21,75); % view(21,50);

 % Explain nuclear mask
      
     
     
     