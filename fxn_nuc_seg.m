%%                                  fxn_nuc_seg.m                        %%
% Alistair Boettiger & Jacques Bothma           Version Completed: 01/28/10          
% Levine Lab                                        Last Modified: 03/10/11

%% Overview
%  I -- image
%  FiltSize -- size of region used in image blurring filter.  (30-50)
% 
% sigmaE -- width of Gaussian excitation filter ~ nuc radius (15)
% simgaI -- width of Gaussina inhibition filter ~nuc diameter (25)
% 
% minN -- minimum nucleus size allowed in pixel area (20 - 100) 

%% Update Log
% updated 07/27/10 to use separate Gaussian filters and difference images
% rather than imorph operations
% Fixed bug 02/03/11: cents were not being computed from the final bw
% segmentation.  
% Updated 03/10/11: Major modifications, new watersheding, fixed bug with
% im2bw, 


function  [bw,cents] = fxn_nuc_seg(I,minN,sigmaE,sigmaI,FiltSize,imblur)


%%  DoG Filter to ID nuclei

% this section needs FiltSize, FiltStd, and nucT

% Four step prep for segmentation. Step 1: enhance contrast. 
% Step 2: Apply difference of gaussian filter to emphasize gaussian looking objects. 
% Step 3: Automatic threshold calculated using Otsu's method.
% Step 4: watershed to split fused nuclei.


% for save data; 
fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';

FiltSize = round(1.3*sigmaI); 

  H = fspecial('disk',imblur); % Filter Kernel       
     I = imfilter(I,H,0); %Apply Filter
    figure(2); clf; imshow(I); 

   I =  adapthisteq(I); %Step 1: enhances the contrast of the grayscale image   
 %   H = - fspecial('log',FiltSize,FiltStd); % Step 2 : Filter Kernel
  
% Advanced Guassian filtering    
   Ex = fspecial('gaussian',FiltSize,sigmaE);
   Ix = fspecial('gaussian',FiltSize,sigmaI);
   
  % Faster method to apply filter -- use straight Gaussians.  
  outE = imfilter(single(I),Ex,'replicate'); 
  outI = imfilter(single(I),Ix,'replicate'); 
  outims = outE - outI;  
 % figure(3); clf; imshow(outims); 
      
  W = watershed(max(outims(:))-outims);
  outims(W==0) = 0; 
  outims = makeuint(outims,16); 
  figure(4); clf; imshow(outims); 
  
   bw = im2bw(outims,graythresh(outims)); %Step 3 : Automatic threshold calculated using Otsu's method.

  
     
%% Count Nuclei

bw = bwareaopen(bw,minN); 
% figure(3); clf; imshow(bw); 

[labeled, nucs] = bwlabel(bw,8);  % count and label nuclei (8-> diagnols count as touch)
S = regionprops(labeled,'Centroid');
cents = reshape([S.Centroid],2,length(S));

% Plotting
CM=label2rgb(labeled, 'hsv', [1,1,1],'shuffle');
CM = 255-CM;
[h,w] = size(I);
Io = uint8(zeros(h,w,3)); 
Io(:,:,1) = makeuint(I,8) - uint8(195*bw);
Io(:,:,2) = makeuint(I,8) - uint8(195*bw);
Io(:,:,3) = makeuint(I,8) - uint8(195*bw);
Io = imadd(Io,CM);
figure(22); clf; imagesc( Io); 
title([num2str(nucs),' total nuclei']); 

