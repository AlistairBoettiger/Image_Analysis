%%                  view3Ddots.m
% Alistair Boettiger                            Date Begun: 02/14/11
%                                               Last Modified: 03/20/11

%% Updates
% Update 03/something to make 3D isosurfaces 
% Updated 03/20/11 to add movie 


clear all;

% folder = '/Volumes/Data/Lab Data/Raw_Data/02-06-11/';
% fname = 'MP10_22C_sna_y_c'; emb = '01';

folder = '/Volumes/Data/Lab Data/Raw_Data/02-17-11/';
datafolder =  '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/';  
% fname = 'MP05_22C_sna_y'; emb = '02';
% fname = 'MP09_22C_hb_y_d'; emb = '01';
fname = 'MGa2x_LacZ_sna'; emb = '05';  yp1 = 200; yp2 = 300; xp1 = 250; xp2 = 350;
%fname = 'MGa_LacZ'; emb = '02';
% fname = 'MGa2x_LacZ_sna_b'; emb = '09';
% fname = 'MP05_22C_sna_y_c'; emb = '05';   yp1 = 100; yp2 = 400; xp1 = 50; xp2 = 350;
stackfolder ='MGa2x/'; % 'MP05_22C/'; %  'MGa1x/';% 


stackdisks = 0;

   Imax = imread([folder,stackfolder,'max_',fname,'_',emb,'.tif']); 
   figure(10); clf; imagesc(Imax);

   [h,w] = size(Imax);
   
%%
  

 % yp1 = 850; yp2 = 1150; xp1 = 1050; xp2 = 1350; % must be square. 
  % yp1 = 1050; yp2 = 1350; xp1 = 1050; xp2 = 1350; % must be square. 
  %  yp1 = 250; yp2 = 550; xp1 = 700; xp2 = 1000; % must be square.   % chosen for 'MGa2x_LacZ_sna'; emb = '02';
   figure(10); clf; imagesc(Imax(xp1:xp2,yp1:yp2,:));
   
     disp(['Coordinates:  ', num2str(xp1), ' : ', num2str(xp2), ',   ' num2str(yp1), ' : ', num2str(yp2) ] );
   
   [hs,ws] = size(Imax(xp1:xp2,yp1:yp2,1));

%%
% filename = [folder,'/',fname];
% Im = lsm_read_mod([filename,'.mat'],str2double(emb),1.5E4); 

Zs = length(Im);

%% Find dots 
%(since it's a small patch this is faster/easier than loading the data).  

tic
disp('finding dots...');

% filtering parameters;  
    sigmaE = 3;% .25  IMPORTANT
    sigmaI = 4; % 3   IMPORTANT
    min_int  = 0.05;    %  5    ;% .05 % not necessary Fix at Zero
    FiltSize = 30;% 
    min_size = 40;% 
    Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
    Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
  
Zs = 40;




% [hs,ws] = size(Il1); 
%Isect = uint16(zeros(hs,ws,Zs));

Isect1 = zeros(hs,ws,Zs);
Isect2 = zeros(hs,ws,Zs);
Inuc = zeros(hs,ws,Zs);
DotData = cell(1,Zs);    
DotMasks = cell(1,Zs); 
im_folder = cell(1,Zs);

for z = 1:Zs % z = 20        
      % Dot finding for channel 1
          im_folder{z} = [folder,stackfolder,fname,'_',emb,'_z',num2str(z),'.tif'];
          Iin_z = imread(im_folder{z}); 
               
          [DotData{z},DotMasks{z}] = dotfinder(Iin_z(xp1:xp2,yp1:yp2,1),Ex,Ix,min_int,min_size);
          bw1 = logical(DotMasks{z}); 
          bnd1 = imdilate(bw1,strel('disk',2)) -bw1;         
         % figure(2); clf; imshow(bnd1);
          mask = double(2*bw1)+bnd1;   
          mask(mask==0)=NaN; 
          mask(mask==1) = 0;       
          % figure(2); clf; imagesc(mask);   
  
          Isect1(:,:,z) = 255*double(Iin_z(xp1:xp2,yp1:yp2,1))/2^16.*mask;
          % figure(2); clf; imagesc(  Isect1(:,:,z));
         
%          % Dot finding for channel 2            
%           [cent2,bw2] = dotfinder(Iin_z(xp1:xp2,yp1:yp2,2),Ex,Ix,min_int,min_size);
%            bw2 = logical(bw2);
%            bnd2 = imdilate(bw2,strel('disk',2)) -bw2;
%          % figure(2); clf; imshow(bnd2);
%           mask = double(2*bw2)+bnd2;   
%           mask(mask==0)=NaN; 
%           mask(mask==1) = 0;       
%           % figure(2); clf; imagesc(mask);          
%          Isect2(:,:,z) = 255*double(Iin_z(xp1:xp2,yp1:yp2,2))/2^16.*mask;

          
    % Nuclear finding            
    In =  Iin_z(xp1:xp2,yp1:yp2,3); % Im{1,z}{3}( xp1:xp2,yp1:yp2 ); 
    nucnoise = double(im2bw(In,.15)); 
    nucnoise(nucnoise==0) = NaN; 
    Inuc(:,:,z) = 255*double(In)/2^16.*nucnoise;
   % figure(2); clf; imagesc(Inuc(:,:,z))   
    
end


%%
mRNAchn = 1; plotZdata = 1; getpreciseZ = 1; consec_layers = 2; ovlap = 2; 
  dotC = CheckDotUpDown(DotData,DotMasks,im_folder,mRNAchn,hs,ws,plotZdata,getpreciseZ,consec_layers,ovlap);


%%
clear Im Iin_z cent2 bw1 bw2 mask

%% show stacks

if stackdisks == 1
    [X,Y] = meshgrid((1:ws)*50,(1:hs)*50);
    Istack = zeros(hs,ws,Zs);
    first = 1; last = Zs;

    nmax = round(max(Inuc(:)));
    c1max = round(max(Isect1(:)));
    c2max = round(max(Isect2(:)));

    % 
    figure(1);  clf;
    colordef black; set(gca,'color','k'); set(gcf,'color','k');

    % Plot nuclei data
    for z=first:last
        In = Inuc(:,:,z) + c1max + c2max; 
        Z = (Zs - z*ones(hs,ws))*340;
        surf(X,Y,Z,In); hold on;
    end
    shading interp;
    alpha(.45);  % make nuclei transparent 



    figure(1);  hold on;
    for z=first:last
        I1 = Isect1(:,:,z) ;
       Z = (Zs - z*ones(hs,ws))*340;
        surf(X,Y,Z,I1); hold on;
    end
    shading interp;

    figure(1);  hold on;
    for z=first:last
        I2 = Isect2(:,:,z) + c1max+1; 
        Z = (Zs - z*ones(hs,ws))*380;
        surf(X,Y,Z,I2); hold on;
    end
    shading interp;

    C2 = zeros(c2max,3);
    for cc = 1:c2max
        C2(cc,:) = [((cc-1)/c2max)^3,((cc-1)/c2max)^.5,((cc-1)/c2max)^2];
    end

    CN = cool(nmax); 
    C1 = hot(c1max);
    %C1 = flipud(1-hot(c1max));


    colormap([0,0,0;C1;C2;CN]); colorbar; caxis([0,nmax+c1max+c2max]); 
    %view(-109,50); 
    view(40,39); axis on;
    set(gca,'FontSize',12);
    zlim([0,Zs*340]);
    xlabel('nanometers');  ylabel('nanometers'); zlabel('nanometers');

    axis off;

end


%%  Prep 3D-isosurfaces view
load([datafolder,fname,'_slidedata'], 'Data'); 


d =Data{str2double(emb),1}.dotC;  
d1 = d(d(:,2)>xp1 & d(:,2)<xp2 & d(:,1)>yp1 & d(:,1)<yp2,:); 

d =Data{str2double(emb),2}.dotC;  
clear Data; 
d2 = d(d(:,2)>xp1 & d(:,2)<xp2 & d(:,1)>yp1 & d(:,1)<yp2,:); 

%%  3D-isosurfaces view 
stp = 3;

x = (1:stp:ws)*50;
y = (1:stp:hs)*50;
 z = Zs*380-380*(1:Zs);
%z = 380*(1:Zs);

[X,Y,Z] = meshgrid(x,y,z);

data = Inuc(1:stp:end,1:stp:end,1:Zs);
data(isnan(data)) = 0;
data = smooth3(data,'box',5);

figure(4); clf;
 patch(isosurface(Y,X,Z,data,65),'FaceColor','blue','EdgeColor','none');
%isonormals(data,h); % alpha .7;
% patch(isocaps(X,Y,Z,data,50),'FaceColor','interp','EdgeColor','none');
%camlight left;


material dull


points1 = [d1(:,1)*50-yp1*50,d1(:,2)*50-xp1*50,380*Zs-380*d1(:,3)]; 
points2 = [d2(:,1)*50-yp1*50,d2(:,2)*50-xp1*50,380*Zs-380*d2(:,3)]; 

V1 = zeros(size(X));  V2 = V1; 

ind1 = cell(Zs,1);
ind2 = cell(Zs,1);
for zs =1:Zs; 
%    xs = rand(1,20)*ws*50/(max(x))*length(x);
%    ys = rand(1,20)*hs*50/(max(y))*length(y);
    layer = d1(:,3)<zs+1 & d1(:,3)>zs -1; 
    xs = points1(layer,1)/(50*ws)*length(x);
    ys = points1(layer,2)/(50*hs)*length(y); 
    ind1{zs} = floor(xs) + length(y)*(floor(ys)) + length(x)*length(y)*(zs-1);
    
    layer = d2(:,3)<zs+1 & d2(:,3)>zs -1; 
    xs = points2(layer,1)/(50*ws)*length(x);
    ys = points2(layer,2)/(50*hs)*length(y); 
    ind2{zs} = floor(xs) + length(y)*(floor(ys)) + length(x)*length(y)*(zs-1);
    
end
inds1 = cell2mat(ind1);
inds2 = cell2mat(ind2);
V1(inds1) = 2;
V2(inds2) = 2;

figure(4); hold on;
patch(isosurface(X,Y,Z,V1,1),'FaceColor','red','EdgeColor','none');
 patch(isosurface(X,Y,Z,V2,1),'FaceColor','green','EdgeColor','none');
camlight left;
light('Position',[1000,1.5E4 8000],'Style','local');
 lighting phong; % view(112,34); 
view(194,80);
%  



% scatter3(points1(:,1),points1(:,2),points1(:,3),'r.');

set(gcf,'color','k'); axis off; 
%% Make movie.
view(0,90);

clear data Imax Isect1 Isect2 nucnoise d X Y Z V1 V2 

 figure(4); 
    oh=findobj(gcf,'type','patch');
for i = 1:36
    if i<24
        rotate(oh,[1 0 0],-4);
    else
        rotate(oh,[0 1 0],10);
    end
    M(i) = getframe(gcf);
end
   

% figure(11); clf;
%      axis off;
%      movie(M,1);
%      
 %    save([datafolder,fname,'_movie'],'M'); 
%%
% load([datafolder,fname,'_movie'],'M'); 
% figure(1); clf;
% movie(M,1);

