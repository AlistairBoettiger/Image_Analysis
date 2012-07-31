%%                                       fxnclean_v2.m
% 
%  Alistair Boettiger                                  Date Begun: 11/25/08
% Levine Lab                                   Version Functional: 03/09/10
% function script.                                  Last Modified: 06/19/10
%
%%  Descrition
% takes an image in 1 to 3 layers, selects largest object, removes
% background objects, and orients largest object major axis to x-axis
%
%% Input variables
%  I          image file 
%  emb_chn    channel from which to get object boundaries. =1 for single
%             channel data
%  objT       Threshold (0 to 1) for detecting objects
%  objC       Optional object imclose strel.  default = 0; 
%  minO       Minimum object size
%
%%  Output variables
% I3         image file, cleaned, oriented, and cropped
%
%% Called by programs:
% imviewer
%
% 
%% Uses codes:
% 
% 
%
%% Updates / Notes
% Updated from fxnclean.m on 6/13/10 to use fast_rotate (running C code)
% instead of imrotate -- inbuilt matlab script.

% Improved speed 06/19/10.  Removed minO, redundant with grab largest
% object.  
% 

%% Active script


 function I3 = fxnclean_v2(I,emb_chn,objT,sigma,imdil,minO)    
%            clear all;
%  load test; 
% I = handles.Io;

FiltSize = 50; 
      % locate object   
   %    figure(2); imshow(I); 
      
   if objT == 0 %loop that allows user to select own threhsold or apply automatic one.
           H = fspecial('gaussian',FiltSize,sigma); % Filter Kernel
            outims = imfilter((I(:,:,emb_chn)),H,'replicate'); %Apply Filter
    %     %  figure(3); clf; imshow(outims)
           %outims=adapthisteq(outims);
           outims=imadjust(outims);   % adapthisteq is really slow
          %  figure(3); clf; imshow(outims)
           bw=im2bw(outims,graythresh(outims));
%           figure(3); clf; imshow(outims)
%          figure(4); clf; imshow(bw);
    else
        bw = im2bw(mat2gray(I(:,:,emb_chn)),objT); %Apply user chosen threshold
    end
       
        bw = imfill(bw,'holes'); % make objects solid   
         %  slow and generally not necessary.  
       
      %   bw = bwareaopen(bw,minO);  % remove objects of less than "minO" pixels   
        bw = imdilate(bw,strel('disk',imdil)); 
        
              figure(22); clf;  % troubleshoot locating of object
              subplot(1,2,1); imshow(I);  
              subplot(1,2,2); imshow(bw); title('embryo selection');      
   
      % choose largest object, remove minors, then orient
         L = bwlabel(bw,8);   % create label matrix of objects in image
         imdata = regionprops(L,'Area', 'Orientation'); % measure properties of objects
         keep = find([imdata.Area]==max([imdata.Area])); % keep largest object
         theta = [imdata.Orientation]; % orientation of objects
         [w, l] = size(L);
         emb = zeros(w,l,3);
         obj_loc = ismember(L,keep); 
         for k=1:3
            emb(:,:,k) = obj_loc; 
         end
         
            
        if strcmp(class(I(:,:,emb_chn)),'uint8')      
         I2 = uint8(double(I).*double(emb)); % deletes all other objects by mult zeros
            % I2=uint8(bsxfun(@times,double(I),outims));       
        elseif strcmp(class(I(:,:,emb_chn)),'uint16')               
           % I2=uint16(bsxfun(@times,double(I),outims));
            I2 = uint16(double(I).*double(emb)); % deletes all other objects by mult zeros 
           else
            disp('Error! unknown input format!')
        end
        % load test2;
         
     %   figure(3); clf; imshow(I2);
   %  save test2;
        
        I2 =  padarray(I2,round([.3*w,.3*l]),0); 
        obj_loc = padarray(obj_loc,round([.3*w,.3*l]),0);       
         I2 = imrotate(I2,360-theta(keep)); % rotate object    
         bw = imrotate(uint8(255*obj_loc),360-theta(keep)); 
         bw = logical(bw);
         
         
         
         bw = imdilate(bw,strel('disk',20)); % add a little buffer;  (may be uncessary, somewhat slow). 
         L = bwlabel(bw,8); 
         %   figure(2); clf; imshow(bw)
         imdata = regionprops(L,'BoundingBox'); % compute boundaries
         c = round([imdata.BoundingBox]);
         try
         I3 = I2(c(2):c(4)+c(2),c(1):c(3)+c(1),:);    
         catch
             I3 = I2(c(2):c(4)+c(2)-1,c(1):c(3)+c(1)-1,:);
         end
         
         figure(20); clf; imshow(I3); title('cleaned,oriented embryo');