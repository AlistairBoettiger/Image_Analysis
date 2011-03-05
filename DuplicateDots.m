%%                              DuplicateDots
% Alistair Boettiger                                   Date Begun: 01/30/11
% Levine Lab                                        Last Modified: 02/01/11


%% Description
% Find replicate dots that are present in the above layer.



%% Updates
% rewritten from scratch 02/01/11 to use direct neighbor finding at full
% resolution.  Reduced resolution to create overlap (original method) is
% not precise enough for my purposes.  

% inds_out: nearest pixels corresponding to centroids
% D2u: centroids of dots unique to layer 2.  

function [inds_out, D2u] = DuplicateDots(ovlap,D1,D2,h,w,plotdata)
fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';

     if isempty(D1)
         inds2 = floor(D2(:,2))+floor(D2(:,1))*h; 
         D2u = D2;
         inds_out = inds2;
     else
            
   


         % convert to indices. 
         inds1 = floor(D1(:,2))+floor(D1(:,1))*h; 
         inds2 = floor(D2(:,2))+floor(D2(:,1))*h; 
         
         inds1(inds1>w*h) = w*h; 
         inds2(inds2>w*h) = w*h; 
         
         R1 = false(h,w); R1(inds1) = 1;
         R2 = false(h,w); R2(inds2) = 1;
         R1z = imdilate(R1,strel('disk',ovlap)); % This is going to slow us down alot  

         %  save([fdata,'/','test2']);
         % load([fdata,'/','test2']);
         
         R2u = R2 - R1z;  R2u(R2u<0) = 0; R2u = logical(R2u); % map of unique dots
        %  R2d = R2 + R1z; R2d(R2d<2) = 0; R2d = logical(R2d); 
        % figure(3); clf; imagesc( R2 + R1z); 
         
         
         R2L = bwlabel(R2u);
      %   R2L = logical(R2u);  % DOESN'T WORK. Foolish matlab autosuggest. 
         R2data = regionprops(R2L,'Centroid'); 
         D2u = reshape([R2data.Centroid],2,length(R2data))'; % vector centroids of unique dots

         % inds_out =  floor(D2u(:,2))+floor(D2u(:,1))*h;  % raster indexed based centroids of unique dots ;
        % inds_out = inds_out'; 

              xx = D2u(:,1);  yy = D2u(:,2);
             inds_out = sub2ind([h,w],xx,yy);
        
        
        %      % Plotting for troubleshooting  
        %         It = uint8(zeros(h,w,3));
        %         It(:,:,1) = uint8(55*R1z);
        %         It(:,:,2) = uint8(155*R2);
        %       figure(1); clf; imshow(It(1:300,1:300,:)); hold on;
        %       plot(D1(:,1),D1(:,2),'m+');
        %       plot(D2u(:,1),D2u(:,2),'co');


        %  % plotting for troubleshooting     
        %      figure(4); clf;  
        %      imshow(Iz(1:300,1:300,:));    hold on;    
        %      plot(D1(:,1),D1(:,2),'go');
        %      plot(D2(:,1),D2(:,2),'y+'); 
    
          
         if plotdata == 1 
             Iz = uint16(zeros(h,w,3));
             Iz(:,:,1) = 5*handles.Im{1,z-1}{handles.mRNAchn1};
             Iz(:,:,3) = 5*handles.Im{1,z}{handles.mRNAchn1};
             figure(5); clf;  
             imshow(Iz(1:300,1:300,:));    hold on;    
             plot(D1(:,1),D1(:,2),'y+');
             plot(D2(:,1),D2(:,2),'co'); 
             plot(D2u(:,1),D2u(:,2),'c.','MarkerSize',10);
         end
    end

