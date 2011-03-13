
%% dotfinder.m 

% Alistair Boettiger                                   Date Begun: 01/21/11
% Levine Lab                                        Last Modified: 03/12/11


%% 
% Adapted from count_all_dots (which is adapted from im_nucdots_exon.

function [cent,bw2,labeled] = dotfinder(I,Ex,Ix,min_int,min_size)
%%     
 % figure(1); clf; imagesc(I);
 %  I = Iin_z( :,:,mRNAchn );%  Im{1,z}{mRNAchn}( xp1:xp2,yp1:yp2 );
 
   bw = im2bw(I,min_int);  % if we're gonna do this, do it first.   
   I(bw==0) =0;
  % figure(2); clf; imagesc(I);

   % Faster method to apply filter -- use straight Gaussians. 
  outE = imfilter(single(I),Ex,'replicate'); 
  outI = imfilter(single(I),Ix,'replicate'); 
  outims = outE - outI;
  
%   figure(10); clf; subplot(2,1,1); imagesc(outE); shading flat;
%   subplot(2,1,2); imagesc(outI); shading flat;
 % figure(11); clf; imagesc(outims); shading flat;  colorbar; 
   


  Iin = makeuint(outims,16); 
  %  figure(11); clf; imagesc(Iin); shading flat;  colorbar; 
  
  
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
        [L,count(t)] = bwlabel(bw); 
    end
  %  figure(1); clf; plot(thresh,count);
    [jnk,im] = max(count);
    imthresh = thresh(im);


%       bw2 = im2bw(Iin,graythresh(Iin));  %  figure(10); clf; imagesc(bw2); 
 
      bw2 = im2bw(Iin,imthresh);
   
   % figure(10); clf; imagesc(bw2); shading flat;colormap hot;
  IO = makeuint(outims.*single(bw2),16); 
%  figure(10); clf; imagesc(IO); shading flat;colormap hot;

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
  
 % figure(10); clf; imagesc(IO); shading flat; colormap hot;
 M = max(IO(:)); 
 L = watershed(double(M-IO));  %  figure(2); clf; imagesc(L); shading flat;
 IO(L==0) = 0;
 %Iw = IO; Iw(L==0)=0; 
 
 % figure(2); clf; imagesc(Iw); shading flat; colormap hot;

% %  OLD WaterShedding;  
%   % figure(2); clf; imshow(bw2);
%   D = -bwdist(~bw2);  
%   L = watershed(D);
%   BW = bw2; BW(L==0)=0; 
%    figure(2); clf; imshow(BW);
  

  bw2 = logical(IO); % do threholding first;
  
  
  % Old threshold last;  
  % bw2 = bw3 & bw2; % Must be above threshold and shape selected by LALI
  % filter and watershedding
  
  
   bw2 = bwareaopen(bw2,min_size);% remove objects less than n pixels in size 
  % figure(10); clf; imshow(bw2);       
  % save([fdata,'/','test2']);
  
% mRNA transcript locating counting
       labeled = bwlabel(bw2,8); % count and label RNAs (8-> diagnols count as touch)   
% labeled =uint16(labeled); 
       R1 = regionprops(labeled,IO,'WeightedCentroid'); % compute mRNA centroids      
       cent = reshape([R1.WeightedCentroid],2,length(R1))';
       
      %  R1 = regionprops(labeled,IO,'Area');
      %  figure(1); clf; hist([R1.Area],100); colordef white;  
      %  figure(11); clf; imagesc(labeled);    length([R1.Area])

end
   
