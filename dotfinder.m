
%% dotfinder.m 

% Alistair Boettiger                            Date Begun: 01/21/11
% Levine Lab 


%% 
% Adapted from count_all_dots (which is adapted from im_nucdots_exon.

function [cent,bw2] = dotfinder(I,alphaE,alphaI,Ex,Ix,min_int,min_size)


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
       
  % Faster method to apply filter -- use straight Gaussians.  
  outE = imfilter(single(I),Ex,'replicate'); 
  outI = imfilter(single(I),Ix,'replicate'); 
  outims = alphaE.*outE - alphaI*outI;
  
  bw2 = im2bw(outims,graythresh(outims)); %Automatic threshold calculated using Otsu's method.  
  % figure(2); clf; imshow(bw2);
  D = -bwdist(~bw2);
  L = watershed(D);
  BW = bw2; BW(L==0)=0; 
  % figure(2); clf; imshow(BW(1:300,1:300));
  
  bw3 = im2bw(I,min_int);   
  bw2 = bw3 & BW; % Must be above threshold and shape selected by LALI filter and watershedding
  bw2 = bwareaopen(bw2,min_size);% remove objects less than n pixels in size 
  % figure(1); clf; imshow(bw2);       
  % save([fdata,'/','test2']);
  % figure(1); clf; imshow(uint8(1-Cell_bnd)); 
  
% mRNA transcript locating counting
       %labeled = bwlabel(bw2,8); % count and label RNAs (8-> diagnols count as touch)    
       labeled = logical(bw2,8); % count and label RNAs (8-> diagnols count as touch)   % logical is faster than bw label.  
        R1 = regionprops(labeled,'Centroid'); % compute mRNA centroids
       cent = reshape([R1.Centroid],2,length(R1))';
% % % plotting results       
%        figure(3); clf; imshow(I);  
%        hold on; scatter(centRNA1(:,1),centRNA1(:,2),'bo','SizeData',90);
%        title(['RNAs = ',num2str(RNAs)],'color','k'); set(gcf,'color','w');      
%        
%          figure(3); clf; imagesc(I); colormap hsv;  colorbar
%        hold on; scatter(centRNA1(:,1),centRNA1(:,2),'bo','SizeData',90);
%        title(['RNAs = ',num2str(RNAs)],'color','k'); set(gcf,'color','w');
     

end
   
