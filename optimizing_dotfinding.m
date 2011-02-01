

 handles.fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';
load([handles.fdata,'/','test']);

%%
  disp('running step 4...'); tic
    alphaE = 1;%  str2double(get(handles.in1,'String')); % alphaE = .955; %
    sigmaE = 2;% str2double(get(handles.in2,'String')); %    sigmaE = 2; %
    alphaI = 1.05;%str2double(get(handles.in3,'String')); %   alphaI = .98; % 
    min_int  = .03;% str2double(get(handles.in4,'String')); %   min_int  = .07; % 
    FiltSize = 25;% str2double(get(handles.in5,'String')); %   FiltSize = 20;% 
    min_size = 10;% str2double(get(handles.in6,'String')); %  min_size = 15; % 
    sigmaI = 3;
   
    % Build the Gaussian Filter   
    Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
    Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
   

    
i = 5;

     I1 = handles.Im{1,i}{handles.mRNAchn1};
     D1 = dotfinder(I1,alphaE,alphaI,Ex,Ix,min_int,min_size);

     I2 = handles.Im{1,i+1}{handles.mRNAchn1};
     D2 = dotfinder(I2,alphaE,alphaI,Ex,Ix,min_int,min_size);
        
    
     
     [h,w] = size(I2);
     
     Iz = uint16(zeros(h,w,3));
     Iz(:,:,1) = 5*I1;
     Iz(:,:,3) = 5*I2;
     
     figure(4); clf;  
     imshow(Iz(1:300,1:300,:));    hold on;    
     plot(D1(:,1),D1(:,2),'go');
     plot(D2(:,1),D2(:,2),'y+'); 
     
     
     
 
     % convert to indices. 
     inds1 = floor(D1(:,2))+floor(D1(:,1))*h; 
     inds2 = floor(D2(:,2))+floor(D2(:,1))*h; 
     R1 = false(h,w); R1(inds1) = 1;
     R2 = false(h,w); R2(inds2) = 1;
     R1z = imdilate(R1,strel('disk',2)); % This is going to slow us down alot  
     
     R2u = R2 - R1z;  R2u(R2u<0) = 0; R2u = logical(R2u); % map of unique dots
     R2L = bwlabel(R2u);
     R2data = regionprops(R2L,'Centroid'); 
     D2u = reshape([R2data.Centroid],2,length(R2data))'; % vector centroids of unique dots
     inds2u =  floor(D2u(:,2))+floor(D2u(:,1))*h;  % raster indexed based centroids of unique dots ;

     
%  %   % Old code: Overlap
%    %  overlap = R2z & R1z;
%   %   R2_unique = bwlabel(overlap); 
%   %   overlap_data = regionprops(R2_unique,'Centroid');  
%   %  R2u = reshape([overlap_data.Centroid],2,length(overlap_data))';

     
      It = uint8(zeros(h,w,3));
      It(:,:,1) = uint8(55*R1z);
      It(:,:,2) = uint8(155*R2);
      figure(1); clf; imshow(It(1:300,1:300,:)); hold on;
      plot(D1(:,1),D1(:,2),'m+');
      plot(D2u(:,1),D2u(:,2),'co');
     
      
                 figure(5); clf;  
     imshow(Iz(1:300,1:300,:));    hold on;    
     plot(D1(:,1),D1(:,2),'y+');
     plot(D2(:,1),D2(:,2),'co'); 
     plot(D2u(:,1),D2u(:,2),'c.','MarkerSize',10);
      
      %%
       
%      % Confirm that dots are right;
%      R = double(R1) + 2*double(R2);  
%      Iz(:,:,2) = uint16(2^16/3*R);  figure(1); clf; imshow(Iz(1:300,1:300,:)); 
   
     toc
     
     scale = 5; 
     inds = DuplicateDots(D1,D2,scale,h,w)'; % custom fxn. 
     
                h2 = ceil(h/scale); w2 = ceil(w/scale); 
                
%                 [ys,xs] = ind2sub([h2,w2],inds_Z{z});
%                xs = floor(xs*h/h2 );
%                ys = floor(ys*h/h2 );
%                hr_ind = sub2ind([h,w],ys,xs); % convert to linear index
%                I1 = false(h,w);
%                I1(hr_ind) = 1; % place all dots on array
               
               
               I2u = false(h2,w2); 
               I2u(inds) = 1;
               I2u = imresize(I2u,h/h2);
               
               Iz(:,:,3) = Iz(:,:,3) + uint16(2^16*I2u);
               
              figure(5); clf;  
     imshow(Iz(1:300,1:300,:));    hold on;    
     plot(D1(:,1),D1(:,2),'y+');
     plot(D2(:,1),D2(:,2),'co'); 
     
     