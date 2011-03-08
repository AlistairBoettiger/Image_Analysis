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

% clear all;
% folder = '/Volumes/Data/Lab Data/Raw_Data/02-17-11/';
% fname = 'MP05_22C_sna_y'; emb = '02';
% %fname = 'MP09_22C_hb_y_d'; emb = '01';
% 
% filename = [folder,'/',fname];
% Im_data = lsm_read_mod([filename,'.mat'],str2double(emb),1.5E4); 
% 
% 
% handles.mRNAchn1 = 1;
% handles.Im = Im_data;

% fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';

% save([fdata,'opt_data'],'handles');
% load([fdata,'opt_data']); 


%%

% zoom in on specific region
 m = .95; % .85;
Zs = length(handles.Im);
[h,w] = size(handles.Im{1,1}{1}); 

xp1= floor(h/2*m); xp2 = floor(h/2*(2-m));
yp1 = floor(w/2*m); yp2 = floor(w/2*(2-m));


    sigmaE = 3;%  IMPORTANT
    sigmaI = 4; % IMPORTANT
    min_int  = 0.02;% .05 % not necessary Fix at Zero
    FiltSize = 30;% 
    min_size = 20;% 
    
   
    % Build the Gaussian Filter   
    Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
    Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
   
Filt = Ex -Ix;
figure(1); clf; imagesc(Filt); colorbar;   

%% Quick view around layer z
  disp('running step 4...'); tic
z = 20;

     I0 = handles.Im{1,z-1}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D0 = dotfinder(I0,Ex,Ix,min_int,min_size);

     I1 = handles.Im{1,z}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D1 = dotfinder(I1,Ex,Ix,min_int,min_size);

     I2 = handles.Im{1,z+1}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D2 = dotfinder(I2,Ex,Ix,min_int,min_size);
     
     I3 = handles.Im{1,z+2}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
     D3 = dotfinder(I3,Ex,Ix,min_int,min_size);
     
    
     [h,w] = size(I2);    
%  plotting    
     Iz = uint16(zeros(h,w,3));
     Iz(:,:,1) = 3*I1 +   I3  + .5*I0;
     Iz(:,:,2) = 1*I2 + 2*I3  +  2*I0;
     Iz(:,:,3) = 3*I2 +   I3  + .4*I0;
  figure(6); clf;  
     imagesc(Iz);    hold on;    
     plot(D1(:,1),D1(:,2),'y+');
     plot(D2(:,1),D2(:,2),'co'); 
     plot(D3(:,1),D3(:,2),'+','color',[1,1,1]);
     plot(D0(:,1),D0(:,2),'+','color',[1,.2,.2]); 

     clear D0 D1 D2 D3  I1 I2 I3 I0;
     
     toc
     %% Find all dots    

% % Quick sample view:
% z = 20;
% Il0 = handles.Im{1,z-1}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
% Il1 = handles.Im{1,z}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
% Il2 = handles.Im{1,z+1}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );
% figure(1); clf; imshow(Il0);


[hs,ws] = size(Iz(:,:,1)); 

% NEEDS to be Double 
Isect1 = zeros(hs,ws,Zs);
% Isect2 = zeros(hs,ws,Zs);
% Inuc = zeros(hs,ws,Zs);



plotdata = 0 ;% don't show 

 DotData = cell(1,Zs);    
  DotMasks = cell(1,Zs); 
 Alldots = uint16(zeros(hs,ws,1)); 

 tic; disp('finding dots...'); 
for z = 1:Zs % z = 20        
      % Dot finding for channel 1
          I1 =  handles.Im{1,z}{handles.mRNAchn1}( xp1:xp2,yp1:yp2 );       
          [cent1,bw1,dL] = dotfinder(I1,Ex,Ix,min_int,min_size);
          % figure(2); clf; imagesc(  I1);
          
          DotData{z} = cent1;
          DotMasks{z} = dL; 
          Alldots(:,:,z) = makeuint(I1,16); % 2 
          % figure(2); clf; imagesc(Alldots(:,:,z)); 
              
          bnd1 = imdilate(bw1,strel('disk',2)) -bw1;         
         % figure(2); clf; imshow(bnd2);
          mask = double(2*bw1)+bnd1;   
          mask(mask==0)=NaN; 
          mask(mask==1) = 0; 
          Isect1(:,:,z) = double(I1).*double(mask);
         %  figure(2); clf; imagesc(mask);   
          % figure(2); clf; imagesc(Isect1(:,:,z) );   colormap hot;
end
toc;
%%
NewDotC = CheckDotUpDown(DotData,DotMasks,Alldots,hs,ws,plotdata);

%% Project all layers

first = 1; last = Zs;
depth_code = jet(last-first+1);   
Alldots_proj = max(Alldots(:,:,first:last),[],3); % perform max project

figure(2); clf; colormap hot; imagesc(Alldots_proj); hold on;
plot(  NewDotC(:,1),NewDotC(:,2),'bo','MarkerSize',10 );

%% 3-D visualization


[X,Y] = meshgrid((1:ws)*50,(1:hs)*50);
Istack = zeros(hs,ws,Zs);
%nmax = round(max(Inuc(:)));
c1max = round(max(Isect1(:)));
%c2max = round(max(Isect2(:)));
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

   %plot3(D2u_Z{z}(:,1)*50,D2u_Z{z}(:,2)*50,((Zs-z)*340)*ones(1,length(D2u_Z{z})),'o','MarkerSize',10,'Color',depth_code(z,:)   ); hold on;
end
plot3(NewDotC(:,1)*50,NewDotC(:,2)*50,340*(Zs-NewDotC(:,3)) ,'.','MarkerSize',20); 

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


C1 = hot(c1max);
%C1 = flipud(1-hot(c1max));
%CN = cool(nmax); 

% colormap([0,0,0;C1;C2;CN]); colorbar; caxis([0,nmax+c1max+c2max]); 
colormap(C1); colorbar; caxis([0,c1max]);
%view(-109,50); 
view(40,-30); axis on;
set(gca,'FontSize',12);
zlim([(Zs-last)*340,(Zs-first+1)*340]);
xlabel('nanometers');  ylabel('nanometers'); zlabel('nanometers');

 
%%  Can we extract the intensities inside the sphere and rplot just the
% sphere with the intensity of the original dot?  

