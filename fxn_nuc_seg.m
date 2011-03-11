%%                                  fxn_nuc_seg.m                        %%
% Alistair Boettiger & Jacques Bothma           Version Completed: 01/28/10          
% Levine Lab                                        Last Modified: 07/27/10

%% Overview
%  I -- image
%  FiltSize -- size of region used in image blurring filter.  (30-50)
% FiltStr -- strength of activation filter (.95 - 1.05)
%  sigmaE -- width of Gaussian excitation filter ~ nuc radius (15)
% simgaI -- width of Gaussina inhibition filter ~nuc diameter (25)
% PercP -- percintile of nuclei considered fused (75 - 100%) 
% minN -- minimum nucleus size allowed in pixel area (20 - 100) 

%% Update Log
% updated 07/27/10 to use separate Gaussian filters and difference images
% rather than imorph operations
% Fixed bug 02/03/11: cents were not being computed from the final bw
% segmentation.  
% 


function  [bw,cents] = fxn_nuc_seg(I,minN,FiltStr,sigmaE,sigmaI,AspectRatio,dilp)


% % %  % Trouble shooting defaults.  
% clear all; load test;  
% FiltSize = 30; % str2double(get(handles.in1,'String'));  %  
% sigmaE = 12; 
% sigmaI = 15;
% a = .96; 
%      

%%  DoG Filter to ID nuclei

% this section needs FiltSize, FiltStd, and nucT

%Three step prep for segmentation. Step 1: enhance contrast. Step 2: Apply
%lorentzian of gaussian filter to emphasize gaussian looking objects. Step
%3: Automatic threshold calculated using Otsu's method.

% figure(3); clf; imshow(I);

% tic
% for save data; 
fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';

FiltSize = round(1.3*sigmaI); 
nucT = 0; % automatic threshold is manditory
% dilp = 3; % dilation/erosion parameter is not modifiable

   I =  adapthisteq(I); %Step 1: enhances the contrast of the grayscale image   
 %   H = - fspecial('log',FiltSize,FiltStd); % Step 2 : Filter Kernel
  
% Advanced Guassian filtering    
  % a = FiltStr; 
   Ex = fspecial('gaussian',FiltSize,sigmaE);
   Ix = fspecial('gaussian',FiltSize,sigmaI);
% H = a.*Ex-Ix;  
%  outims = imfilter(double(I),H,0); %Apply Filter       
%    figure(2); clf; surf(H); shading flat; camlight left; lighting gouraud;
   
  % Faster method to apply filter -- use straight Gaussians.  
  outE = imfilter(single(I),Ex,'replicate'); 
  outI = imfilter(single(I),Ix,'replicate'); 
  outims = outE - outI;  

 % figure(3); clf; imshow(outims); 
      
  W = watershed(max(outims(:))-outims);
  outims(W==0) = 0; 
  outims = makeuint(outims,16); 
  figure(4); clf; imshow(outims); 
  
   % Set negative values to zero and scale to span 16 bit color depth
   % outims(outims<0) = 0; 
   % outims=uint16(outims/max(outims(:))*(2^16-1));

    if nucT == 0 %loop that allows user to select own threhsold or apply automatic one.
         bw = im2bw(outims,graythresh(outims)); %Step 3 : Automatic threshold calculated using Otsu's method.
    else
        bw = im2bw(outims,nucT); %Apply user chosen threshold
    end
  
    
%       D = -bwdist(~bw);
%       L = watershed(D);
%       bw(L==0)=0; 

    
    
% % Plotting stuff    
    % handles.bw = bw;  % Output binary image to next step
  %  L = bwlabeln(bw,8); % Label the unique regions in the thresholded image
    % CM=label2rgb(L, 'jet', [1,1,1],'shuffle');
    % DI = uint16(bsxfun(@times,double(CM)/255,double(I)));
%   figure(3), imshow(DI,'Border','tight','InitialMagnification',100); % maxwindow
% toc
%% Count Nuclei
%  this section needs PercP, dilp, and minN

% nuc_bw = bw;
% 
% %  clean up large and fused nuclei
% L = bwlabeln(nuc_bw,8); % label
% S = regionprops(L,'Perimeter','Area','Centroid','MajorAxisLength','MinorAxisLength'); % measure areas
% bwjn = ismember(L,find([S.MajorAxisLength] < AspectRatio*[S.MinorAxisLength] )); % map of all unjoined nuclei
% bwj =  ismember(L,find([S.MajorAxisLength] > AspectRatio*[S.MinorAxisLength] ));% map of all joined nuclei
% bwj = imerode(bwj,strel('disk',dilp)); %errode joined nuclei
% bwj = imdilate(bwj,strel('disk',dilp-1));
% 

% 
% 
% % save([fdata,'/','test2']); % figure(3); clf; imshow(bw);
% 
% bw = logical(bwj + bwjn); 
bw = bwareaopen(bw,minN); 
% figure(3); clf; imshow(bw); 

[labeled, nucs] = bwlabel(bw,8);  % count and label nuclei (8-> diagnols count as touch)

S = regionprops(labeled,'Centroid');
cents = reshape([S.Centroid],2,length(S));
CM=label2rgb(labeled, 'hsv', [1,1,1],'shuffle');
CM = 255-CM;

 [h,w] = size(I);
Io = uint8(zeros(h,w,3)); 
Io(:,:,1) = makeuint(I,8) - uint8(195*bw);
Io(:,:,2) = makeuint(I,8) - uint8(195*bw);
Io(:,:,3) = makeuint(I,8) - uint8(195*bw);
Io = imadd(Io,CM);
%figure(21); clf; imagesc( Io);




% Io(:,:,1) = .5*(makeuint(I,8)-CM(:,:,1))  +.5*makeuint(I,8).*CM(:,:,1);
% Io(:,:,2) = .5*(makeuint(I,8)-CM(:,:,2))  +.5*makeuint(I,8).*CM(:,:,2);
% Io(:,:,3) = .5*(makeuint(I,8)-CM(:,:,3))  +.5*makeuint(I,8).*CM(:,:,3);


% Io(:,:,1) =2*makeuint(I,8) - .5*CM(:,:,1);
% Io(:,:,2) =2*makeuint(I,8) - .5*CM(:,:,2);
% Io(:,:,3) =2*makeuint(I,8) - .5*CM(:,:,3);
figure(22); clf; imagesc( Io);  title([num2str(nucs),' total nuclei']); 

% toc 
% save test2;

