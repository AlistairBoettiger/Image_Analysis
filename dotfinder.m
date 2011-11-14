
%% dotfinder.m 

% Alistair Boettiger                                   Date Begun: 01/21/11
% Levine Lab                                        Last Modified: 06/06/11


%% 
% Adapted from count_all_dots (which is adapted from im_nucdots_exon)
%% Modified:
% modified 3/12/11 to improve method
% modified 03/13/11 to reuse more variable names more and reduce creating
% mulitple large arrays. 
% modified 06/06/11 to output dot_labels from the label matrix (in place of
% the whole label matrix), the index of the center pixels (in addition to
% the pixels themselves) and the intensity values at the center pixels.

function [dot_labels,cent,inds,ints] = dotfinder(I,Ex,Ix,min_int,min_size)
%%     
 % figure(1); clf; imagesc(I);
 %  I = Iin_z( :,:,mRNAchn );    %  Im{1,z}{mRNAchn}( xp1:xp2,yp1:yp2 );
%  
%  xp1 = 640; yp1 = 640; xp2 = 840; yp2 = 840; 
%  xp1 = 340; yp1 = 340; xp2 = 540; yp2 = 540; 
%   z = 20;  Iin_z = imread(im_folder{z}); I = Iin_z(xp1:xp2,yp1:yp2,1);
% figure(9); clf; imagesc(3*I); ; colormap hot; set(gcf,'color','k');  colorbar; 
  
  
   bw = im2bw(I,min_int);  % if we're gonna do this, do it first.   
   Iin = I;
   Iin(bw==0) =0;
   clear bw; 
  % figure(2); clf; imagesc(Iin); colormap hot; set(gcf,'color','k');  colorbar; 
    % figure(8); clf; imagesc(bw); colormap hot; 
   % Faster method to apply filter -- use straight Gaussians. 
  outE = imfilter(single(Iin),Ex,'replicate'); 
  outI = imfilter(single(Iin),Ix,'replicate'); 
  Iin = makeuint(outE - outI,16);
  
%   figure(10); clf; subplot(2,1,1); imagesc(outE); shading flat;
%   subplot(2,1,2); imagesc(outI); shading flat;
 % figure(8); clf; imagesc(1.2*Iin); shading flat;  colorbar; colormap hsv; set(gcf,'color','k');     
    
  
  % Max dots selection of threshold.  
  % Makes 2048 x 2048 twice as long (700+ seconds compared to 350s). 
    % Use only a subdomain to pick the threshold:
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
  %  figure(1); clf; plot(thresh,count,'linewidth',2); set(gcf,'color','k'); xlabel('theta'); ylabel('Number of objects');  
    [jnk,im] = max(count);
    imthresh = thresh(im);
    bw2 = im2bw(Iin,imthresh);

%       bw2 = im2bw(Iin,graythresh(Iin));
      
   Iin = Iin.*uint16(bw2);  
% figure(10); clf; imagesc(bw2); shading flat; colormap hot; set(gcf,'color','k');  
%  figure(10); clf; imagesc(Iin); shading flat;colormap hot;

   
   %Iin = makeuint(outims.*single(bw2),16); 
  

% ----------------------------------------------------------------- %
% % BUG BAIT!! 
%  %  if given an input of type single or double im2bw thresholds at x>T
%  is a 1 and x <T  is a zero, regardless of the range of x. T is however
%  restricted to be between 0 and 1.  
%   bw_old = im2bw(outims,graythresh(Iin)); 
%   I_old = I.*single(bw_old); figure(11); clf; imagesc(bwO); 
%     figure(11); clf; imagesc(I_old); shading flat; colormap hot; 
%  regionprops must also be given uint16 intensity matrices or doubles and
%  singles in the range [0,1].
% ----------------------------------------------------------------- %  
  
 %figure(10); clf; imagesc(Iin); shading flat; colormap hot; set(gcf,'color','k');    
 M = max(Iin(:)); 
 L = watershed(double(M-Iin));  %  figure(2); clf; imagesc(L); colormap lines; shading flat;
 Iin(L==0) = 0;
 
% %  OLD WaterShedding;  
%   % figure(2); clf; imshow(bw2);
%   D = -bwdist(~bw2);  
%   L = watershed(D);
%   BW = bw2; BW(L==0)=0; 
%    figure(2); clf; imshow(BW);
   
    % Old threshold last;  
  % bw2 = bw3 & bw2; % Must be above threshold and shape selected by LALI
  % filter and watershedding
  
  
   bw2 = logical(Iin); 
   bw2 = bwareaopen(bw2,min_size);% remove objects less than n pixels in size 
  % figure(10); clf; imagesc(bw2); shading flat; colormap hot; set(gcf,'color','k');       
  % save([fdata,'/','test2']);
  
% mRNA transcript locating counting
       labeled = bwlabel(bw2,8); % count and label RNAs (8-> diagnols count as touch)  
    %    figure(2); clf; imagesc(labeled); colormap lines; shading flat;
% labeled =uint16(labeled); 
       R1 = regionprops(labeled,Iin,'WeightedCentroid'); % compute mRNA centroids      
       cent = reshape([R1.WeightedCentroid],2,length(R1))';
       
       
  % % RECORD indices of each dot and intensity at centroid of each dot
       inds = floor(cent(:,2))+floor(cent(:,1))*h;  % indices in this layer  
       inds(inds>w*h) = w*h; 
       ints = [I(inds);0]; 
       dot_labels = labeled(inds); 
       % the leading zero is a dummy which allows pixel values to be set to
       % 0 while referencing members of ints. 
         
         
   
       
%        
%        figure(10); clf; imagesc(uint16(bw2).*Iin); hold on;
%        plot(cent(:,1),cent(:,2),'bo');
%        figure(7); clf; imagesc(3*I); hold on; 
%         shading flat; colormap hot; set(gcf,'color','k'); 
%        plot(cent(:,1),cent(:,2),'bo');
       
      %  R1 = regionprops(labeled,IO,'Area');
      %  figure(1); clf; hist([R1.Area],100); colordef white;  
      %  figure(11); clf; imagesc(labeled);    length([R1.Area])

end
   
