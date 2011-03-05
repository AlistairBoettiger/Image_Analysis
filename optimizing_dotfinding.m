%%                          optimizing _dotfinding.m
% Alistair Boettiger                                   Date Begun: 01/31/11
% Levine Lab 

% Just a testing script to optimize dot finding algorithm and parameters.
% This is a faster testbed then running the whole stack analysis as in
% imsinglemolecule.fig
% 
% 

%  handles.fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';
% load([handles.fdata,'/','test']);

folder = '/Volumes/Data/Lab Data/Raw_Data/02-17-11/';
fname = 'MP05_22C_sna_y'; emb = '02';
%fname = 'MP09_22C_hb_y_d'; emb = '01';

filename = [folder,'/',fname];
% Im_data = lsm_read_mod([filename,'.mat'],str2double(emb),1.5E4); 

%%
handles.mRNAchn1 = 1;
handles.Im = Im_data;
Im = Im_data;



% zoom in on specific region
 m = .9; % .85;
Zs = length(Im_data);
[h,w] = size(Im_data{1,1}{1}); 

xp1= floor(h/2*m); xp2 = floor(h/2*(2-m));
yp1 = floor(w/2*m); yp2 = floor(w/2*(2-m));


    alphaE = 1;% 
    sigmaE = 3;% 
    alphaI = 1.0;%
    min_int  = 0.01;% .05
    FiltSize = 25;% 
    min_size = 10;% 
    sigmaI = 3.2;
   
    % Build the Gaussian Filter   
    Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
    Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
   
Filt = alphaE*Ex - alphaI*Ix;
%figure(1); clf; imagesc(Filt); colorbar;   

%%
  disp('running step 4...'); tic
z = 20;

     I0 = handles.Im{1,z-1}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D0 = dotfinder(I0,alphaE,alphaI,Ex,Ix,min_int,min_size);


     I1 = handles.Im{1,z}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D1 = dotfinder(I1,alphaE,alphaI,Ex,Ix,min_int,min_size);

     I2 = handles.Im{1,z+1}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D2 = dotfinder(I2,alphaE,alphaI,Ex,Ix,min_int,min_size);
     
     I3 = handles.Im{1,z+2}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D3 = dotfinder(I3,alphaE,alphaI,Ex,Ix,min_int,min_size);
     
       I4 = handles.Im{1,z+3}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D4 = dotfinder(I4,alphaE,alphaI,Ex,Ix,min_int,min_size);
    
     
     [h,w] = size(I2);
     
 
     
     
     
     
 
     % convert to indices. 
     inds1 = floor(D1(:,2))+floor(D1(:,1))*h; 
     inds2 = floor(D2(:,2))+floor(D2(:,1))*h; 
     inds1(inds1>w*h) = w*h; 
     inds2(inds2>w*h) = w*h; 
     
     R1 = false(h,w); R1(inds1) = 1;
     R2 = false(h,w); R2(inds2) = 1;
     R1z = imdilate(R1,strel('disk',2)); % This is going to slow us down alot  
     
     R2u = R2 - R1z;  R2u(R2u<0) = 0; R2u = logical(R2u); % map of unique dots
     R2L = bwlabel(R2u);
     R2data = regionprops(R2L,'Centroid'); 
     D2u = reshape([R2data.Centroid],2,length(R2data))'; % vector centroids of unique dots
     inds2u =  floor(D2u(:,2))+floor(D2u(:,1))*h;  % raster indexed based centroids of unique dots ;
 
%  plotting    
     Iz = uint16(zeros(h,w,3));
     Iz(:,:,1) = 5*I1 + 2*I4+    I3;
     Iz(:,:,2) =        3*I4 + 2*I3;
     Iz(:,:,3) = 5*I2+  2*I4+    I3;
  figure(5); clf;  
     imshow(Iz);    hold on;    
     plot(D1(:,1),D1(:,2),'y+');
     plot(D2(:,1),D2(:,2),'co'); 
     plot(D2u(:,1),D2u(:,2),'c.','MarkerSize',10);
       plot(D3(:,1),D3(:,2),'+','color',[1,1,1]);
     plot(D4(:,1),D4(:,2),'+','color',[.5,1,.25]); 
      
   
     toc
     %%
     
     



[h,w] = size(Im{1,1}{1});



% % Quick sample view:
% z = 20;
% Il0 = Im{1,z-1}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
% Il1 = Im{1,z}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
% Il2 = Im{1,z+1}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
% figure(1); clf; imshow(Il0);



[hs,ws] = size(Iz(:,:,1)); 
%Isect = uint16(zeros(hs,ws,Zs));

Isect1 = zeros(hs,ws,Zs);
Isect2 = zeros(hs,ws,Zs);
Inuc = zeros(hs,ws,Zs);
Cents = cell(1,Zs); 

        inds_Z = cell(Zs,1);
        D2u_Z = cell(Zs,1); % store 
        plotdata = 0 ;% don't show 

        
        ovlap = 5; 
        
 DotData = cell(1,Zs);    
  DotMasks = cell(1,Zs); 
 Alldots = uint16(zeros(hs,ws,1)); 

 tic; disp('finding dots...'); 
for z = 1:Zs % z = 20        
      % Dot finding for channel 1
          I1 =  Im{1,z}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
          [cent1,bw1,dL] = dotfinder(I1,alphaE,alphaI,Ex,Ix,min_int,min_size);
         
          % figure(2); clf; imagesc(  Isect1(:,:,z));
          
          DotData{z} = cent1;
          DotMasks{z} = dL; 
           Alldots(:,:,z) = I1; % 2 

 % % Old method: DuplicateDots        
%             if z == 1 % z = 5;
%                 D1 = []; % there is no previous layer if z = 1; 
%             else
%                 D1 = DotData{z-1}; % all dots in the previous layer
%             end
%                 D2 = DotData{z}; % all dots in this layer
%                 [inds_Z{z}, D2u_Z{z}] = DuplicateDots(ovlap,D1,D2,hs,ws,plotdata); % custom fxn. 
  %               % should have updated to hs ws earlier ? ! 
end
toc;
%%
    tic; disp('checking dots...'); 
for z=1:Zs

         % New Method: Check layer above and below for same dot;        
            if z == 1 % z = 5;
                D0 = NaN*ones(1,2); % there is no previous layer if z = 1; 
            else
                D0 = DotData{z-1};  % all dots in the previous layer
            end
            
            D1 =DotData{z};  % all dots in this layer
            
            if z==Zs
                D2 = NaN*ones(1,2);
            else
                D2 = DotData{z+1}; 
            end
            [inds_Z{z}, D2u_Z{z},R3d] = CheckDotUpDown(ovlap,D0,D1,D2,hs,ws,plotdata);
          
          
          bw1 =  imdilate(R3d,strel('disk',5));  
          bnd1 = imdilate(bw1,strel('disk',2)) -bw1;         
         % figure(2); clf; imshow(bnd2);
          mask = double(2*bw1)+bnd1;   
          mask(mask==0)=NaN; 
          mask(mask==1) = 0;       
          % figure(2); clf; imagesc(mask);    
          I1 =  Im{1,z}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
          Isect1(:,:,z) = 255*double(I1)/2^16.*mask;
            
            
                % Returns indices in layer 2 that are not also in layer 1. 
         
         % For projecting all dots into single plane 
            
          
          
          
          
%           % Dot finding for channel 2
%           I2 = Im{1,z}{2}( xp1:xp2,yp1:yp2 );
%           [cent2,bw2] = dotfinder(I2,alphaE,alphaI,Ex,Ix,min_int,min_size);
%           bnd2 = imdilate(bw2,strel('disk',2)) -bw2;         
%          % figure(2); clf; imshow(bnd2);
%           mask = double(2*bw2)+bnd2;   
%           mask(mask==0)=NaN; 
%           mask(mask==1) = 0;       
%           % figure(2); clf; imagesc(mask);          
%           Isect2(:,:,z) = 255*double(I2)/2^16.*mask;
%           
%     % Nuclear finding            
%     In = Im{1,z}{3}( xp1:xp2,yp1:yp2 ); 
%     nucnoise = double(im2bw(In,.15)); 
%     nucnoise(nucnoise==0) = NaN; 
%     Inuc(:,:,z) = 255*double(In)/2^16.*nucnoise;
%     figure(2); clf; imagesc(Inuc(:,:,z));
    
    
end
toc;

%%

first = 1; last = Zs;
depth_code = jet(last-first+1);   
Alldots_proj = max(Alldots(:,:,first:last),[],3); % perform max project


figure(2); clf; colormap hot; imagesc(Alldots_proj); hold on;
for z=first:last
        
     %   plot(  DotData{z}(:,1),DotData{z}(:,2),'.','MarkerSize',10,'Color',depth_code(z,:)); hold on;
        try
         plot(  D2u_Z{z}(:,1),D2u_Z{z}(:,2),'o','MarkerSize',5,'Color',depth_code(z,:)); hold on;
        catch me
            disp(me.message);
        end
   
  %  end
%     figure(11);  hold on;
%      scatter3(  cent1(:,1)*50,cent1(:,2)*50,(Zs-z*340)*ones(1,length(cent1)),'r.','SizeData',100   );
end
%%


[X,Y] = meshgrid((1:ws)*50,(1:hs)*50);
Istack = zeros(hs,ws,Zs);
nmax = round(max(Inuc(:)));
c1max = round(max(Isect1(:)));
c2max = round(max(Isect2(:)));
figure(1);  clf;
colordef black; set(gca,'color','k'); set(gcf,'color','k');

% % Plot nuclei data
% for z=first:last
%     In = Inuc(:,:,z) + c1max + c2max; 
%     Z = (Zs - z*ones(hs,ws))*340;
%     surf(X,Y,Z,In); hold on;
% end
% shading interp;
% alpha(.45);  % make nuclei transparent 

figure(1); clf;  
for z=first:last
    figure(1);  hold on;
    I1 = Isect1(:,:,z) ;
   Z = (Zs - z*ones(hs,ws))*340;
    surf(X,Y,Z,I1); hold on;
   % if z>2 && z<10
       % plot3(  DotData{z}(:,1)*50,DotData{z}(:,2)*50,((Zs-z)*340)*ones(1,length(DotData{z})),'.','MarkerSize',10,'Color',depth_code(z,:)   );
  %  end
   plot3(D2u_Z{z}(:,1)*50,D2u_Z{z}(:,2)*50,((Zs-z)*340)*ones(1,length(D2u_Z{z})),'o','MarkerSize',10,'Color',depth_code(z,:)   ); hold on;
end
shading interp;

% figure(1);  hold on;
% for z=first:last
%     I2 = Isect2(:,:,z) + c1max+1; 
%     Z = (Zs - z*ones(hs,ws))*340;
%     surf(X,Y,Z,I2); hold on;
% end
% shading interp;
% 
% C2 = zeros(c2max,3);
% for cc = 1:c2max
%     C2(cc,:) = [((cc-1)/c2max)^3,((cc-1)/c2max)^.5,((cc-1)/c2max)^2];
% end

CN = cool(nmax); 
C1 = hot(c1max);
%C1 = flipud(1-hot(c1max));


% colormap([0,0,0;C1;C2;CN]); colorbar; caxis([0,nmax+c1max+c2max]); 
colormap(C1); colorbar; caxis([0,c1max]);
%view(-109,50); 
view(40,39); axis on;
set(gca,'FontSize',12);
zlim([(Zs-last)*340,(Zs-first+1)*340]);
xlabel('nanometers');  ylabel('nanometers'); zlabel('nanometers');

 