
%% dotfinder.m 

% Alistair Boettiger                            Date Begun: 01/21/11
% Levine Lab 


%% 
% Adapted from count_all_dots (which is adapted from im_nucdots_exon.

function [cent,bw2,labeled] = dotfinder(I,alphaE,alphaI,Ex,Ix,min_int,min_size)


%  % Optional display       
%    %~~~~~~~~Plot Filter in Upper Right Corner of Image~~~~~~~~~~~%
%        H = alphaE.*Ex-alphaI*Ix;  
%        norm_fact = max(max(H));
%        Hn = 255*H./norm_fact; 
%        sI = Hn;  sI(Hn>0) = 0; sI = uint8(-sI);  
%        sE = uint8(Hn);  sE(Hn<0) = 0; 
% 
%        
%        [h,w] = size(I2); 
%        Xc = round(w/10);
%        Yc = round(h/10);
% 
%        I3 = uint8(zeros(h,w,3));
%        I3(:,:,2) = I2;
%        I3(Yc+1:Yc+FiltSize,Xc+1:Xc+FiltSize,1) = sI;
%        I3(Yc+1:Yc+FiltSize,Xc+1:Xc+FiltSize,3) = sE;  
%      %   figure(2); clf; imshow(I3); 
%   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  
  % outims = imfilter(double(I),H,0); %Apply Filter    ## Old way to apply
  % filter  
       
 
  % I = I1; 
  

   % Faster method to apply filter -- use straight Gaussians. 
  outE = imfilter(single(I),Ex,'replicate'); 
  outI = imfilter(single(I),Ix,'replicate'); 
  outims = alphaE.*outE - alphaI*outI;
  
%   figure(10); clf; subplot(2,1,1); imagesc(outE); shading flat;
%   subplot(2,1,2); imagesc(outI); shading flat;
%   figure(11); clf; imagesc(outims); shading flat;  colorbar; 
   

 % Automatic threshold calculated using Otsu's method.  
  Iin = makeuint(outims,16); 
  bw2 = im2bw(Iin,graythresh(Iin));  %  figure(10); clf; imagesc(bw2); 
  IO = outims.*single(bw2); 
%  figure(10); clf; imagesc(IO); shading flat;colormap hot;

% ----------------------------------------------------------------- %
% % Matlab BUG FOUND !! 
%  %  if given an input of type single or double im2bw thresholds at x>T
%  is a 1 and x <T  is a zero, regardless of the range of x. T is however
%  restricted to be between 0 and 1.  
%   bw_old = im2bw(outims,graythresh(Iin)); 
%   I_old = I.*single(bw_old); figure(11); clf; imagesc(bwO); 
%     figure(11); clf; imagesc(I_old); shading flat; colormap hot; 
% ----------------------------------------------------------------- %  
  
 % figure(10); clf; imagesc(IO); shading flat; colormap hot;
 M = max(IO(:)); 
 L = watershed(M-double(IO));  %  figure(2); clf; imagesc(L); shading flat;
 Iw = IO; Iw(L==0)=0; 
 % figure(2); clf; imagesc(Iw); shading flat; colormap hot;

% %  OLD WaterShedding;  
%   % figure(2); clf; imshow(bw2);
%   D = -bwdist(~bw2);  
%   L = watershed(D);
%   BW = bw2; BW(L==0)=0; 
%    figure(2); clf; imshow(BW);
  

  bw2 = logical(Iw); % do threholding first;
  
  
  % Old threshold last;  
   bw3 = im2bw(I,min_int);   
   bw2 = bw3 & bw2; % Must be above threshold and shape selected by LALI
  % filter and watershedding
  
  
  bw2 = bwareaopen(bw2,min_size);% remove objects less than n pixels in size 
  % figure(1); clf; imshow(bw2);       
  % save([fdata,'/','test2']);
  % figure(1); clf; imshow(uint8(1-Cell_bnd)); 
  
% mRNA transcript locating counting
       labeled = bwlabel(bw2,8); % count and label RNAs (8-> diagnols count as touch)   
       labeled =uint16(labeled); 
       R1 = regionprops(labeled,IO,'WeightedCentroid'); % compute mRNA centroids
      cent = reshape([R1.WeightedCentroid],2,length(R1))';
       
%        % simple centroid
%         R1 = regionprops(labeled,'Centroid'); % compute mRNA centroids
%         cent = reshape([R1.Centroid],2,length(R1))';
       
       
% % % plotting results       
%        figure(3); clf; imshow(I);  
%        hold on; scatter(centRNA1(:,1),centRNA1(:,2),'bo','SizeData',90);
%        title(['RNAs = ',num2str(RNAs)],'color','k'); set(gcf,'color','w');      
%        
%          figure(3); clf; imagesc(I); colormap hsv;  colorbar
%        hold on; scatter(centRNA1(:,1),centRNA1(:,2),'bo','SizeData',90);
%        title(['RNAs = ',num2str(RNAs)],'color','k');
%        set(gcf,'color','w');
end
   
