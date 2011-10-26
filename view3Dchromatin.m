%%                         view3Dchromatin.m
% Alistair Boettiger                                   Date Begun: 02/14/11
%                                                   Last Modified: 06/17/11

%% Updates



clear all;

% folder = '/Volumes/Data/Lab Data/Raw_Data/02-06-11/';
% fname = 'MP10_22C_sna_y_c'; emb = '01';

% folder = '/Volumes/Data/Lab Data/Raw_Data/02-17-11/';
% datafolder =  '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/';  
% stackfolder = 'CadN/'; % 'MP09_22C/'; % 
 % fname ='MP09_22C_hb_y_a'; emb = '01' ;  xp1 = 600; yp1 = 660; xp2 = xp1+100;  yp2 = yp1+100; 
 % fname =  'CadN';  emb = '06';     xp1 = 1275; yp1 = 1050; xp2 = xp1+100;  yp2 = yp1+100;  


datafolder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-11/'; %
folder = '/Volumes/Data/Lab Data/Raw_Data/2011-05-22/'; % /02-17-11/'; %%   %
stackfolder =  's01_MP09/';%  's03_MP02/';%   's02_MP01/';%      's10_bcd1x/';%  's11_bcd6x/'; %'s14_comp_cntrl/'; % 's12_cntrl_2label/'; %'MP02_22C/'; %'MP01_22C/'; % 'MGa1x/'; % 'MP10_22C/'; %'MP05_22C/'; %'YW_ths_sog/'; % 'MP10_22C/'; %  % 'MP09_22C/'; % 'MGa2x/'; % 'MGa1x/'; % 'MGa2x/'; % 'MP10_22C_sna_y_c/'; %
fname =  's01_MP09_Hz_22C';  emb = '09';  %  's03_MP02_Hz_22C'; emb = '01'; %  's02_MP01_Hz_22C'; emb = '04'; %    's10_bcd1x';% 's11_bcd6x'; % 's14_comp_cntrl'; Es =1; % 's12_cntrl_2label'; Es = 1; % 'MP09_22C_hb_y_f'; Es = 7; %  'MP02_22C_hb_y'; Es = 9; % 'MP02_22C_hb_y_b'; Es = 10; %  % 'MP01_22C_hb_y_f'; Es = 12; % 'MP01_22C_hb_y_c'; Es = 10; % 'MP01_22C_hb_y'; Es = 13; % 'MGa1x_LacZ_b'; Es = 12; %  'MP10_22C_sna_y_e'; Es = 12; %  'MP05_22C_sna_y_c'; Es =7; %  'MP10_22C_sna_y_d3'; Es = 1;  %'YW_ths_sog'; Es = 12;  % % 'MP09_22C_hb_y_e'; Es = 10; % 'MP09_22C_hb_y_d'; Es=11; % 'MGa2x_LacZ_sna_b'; Es = 10; % 'MP10_22C_sna_y_d';   % 'MGa_LacZ'; %'MGa2x_LacZ_sna'; %'MP10_22C_sna_y_c'; old_lab = 1;  % 'MP05_22C_sna_y'; old_lab = 1; % 
 
 
 stackdisks = 0;  zspace = 330;
 shownucs = 0;  disp_nascent = 0;

   Imax = imread([folder,stackfolder,'max_',fname,'_',emb,'.tif']); 

   
%%
%xp1 = 600; yp1 = 660; xp2 = xp1+100;  yp2 = yp1+100;   nt = 5; 
% xp1 = 1850; yp1 = 990; xp2 = xp1+100;  yp2 = yp1+100;   nt = 5; 

   xp1 = 145; yp1 = 1660; xp2 = xp1+260;  yp2 = yp1+260; nt = 5; % MP09_09;
  % xp1 = 1245; yp1 = 1260; xp2 = xp1+260;  yp2 = yp1+260;  nt = 15;  % MP01_04 nascent;
  % xp1 = 200; yp1 = 1200; xp2 = xp1+260;  yp2 = yp1+260;   nt = 15; % MP01_04; dots
 % xp1 = 1130; yp1 = 700; xp2 = xp1+260;  yp2 = yp1+260;   nt = .5; % MP02_01; 
   
   % Imax(:,:,2) = 10*Imax(:,:,2);
  
   figure(10); clf; imagesc(Imax); colordef black; set(gcf,'color','k');
   hold on; plot([yp1,yp2,yp2,yp1,yp1],[xp1,xp1,xp2,xp2,xp1],'w');
 
   figure(9); clf; imagesc(Imax(xp1:xp2,yp1:yp2,:)); 
   disp(['Coordinates:  ', num2str(xp1), ' : ', num2str(xp2), ',   ' num2str(yp1), ' : ', num2str(yp2) ] );   
   [hs,ws] = size(Imax(xp1:xp2,yp1:yp2,1));
   [h,w] = size(Imax);




%% Load data into stacks
tic
disp('loading data...');
Zs = 50;

im_folder = cell(Zs,1); 
Isect1 = zeros(hs,ws,Zs);
Isect2 = zeros(hs,ws,Zs);
Inuc = zeros(hs,ws,Zs);

for z = 1:Zs % z = 20     
    try
          im_folder{z} = [folder,stackfolder,fname,'_',emb,'_z',num2str(z),'.tif'];
          Iin_z = imread(im_folder{z});                
    catch err
        disp(err.message); 
        break
    end
          
    % Nuclear finding            
    In =  Iin_z(xp1:xp2,yp1:yp2,3); % Im{1,z}{3}( xp1:xp2,yp1:yp2 ); 
    nucnoise = double(im2bw(In,.15)); 
    nucnoise(nucnoise==0) = NaN; 
    Inuc(:,:,z) = 255*double(In)/2^16.*nucnoise;
   % figure(2); clf; imagesc(Inuc(:,:,z))   
   
   Isect1(:,:,z) =  Iin_z(xp1:xp2,yp1:yp2,1); 
   Isect2(:,:,z) = Iin_z(xp1:xp2,yp1:yp2,2);
    
end

toc

%%  load mRNA data Prep 3D-isosurfaces view
% load([datafolder,fname,'_slidedata'], 'Data'); 
load([datafolder,fname,'_',emb,'_chn1','_data.mat']);

d = dotC;  
d1 = d(d(:,2)>xp1 & d(:,2)<xp2 & d(:,1)>yp1 & d(:,1)<yp2,:); 

load([datafolder,fname,'_',emb,'_chn2','_data.mat']);
d = dotC;  
d2 = d(d(:,2)>xp1 & d(:,2)<xp2 & d(:,1)>yp1 & d(:,1)<yp2,:); 



%%  3D-isosurfaces view, Prep
stp = 1;
Zmax = 50; 

x = (1:stp:ws)*50;
y = (1:stp:hs)*50;
 z = Zs*zspace-zspace*(1:Zmax);

[X,Y,Z] = meshgrid(x,y,z);

data = Inuc(1:stp:end,1:stp:end,1:Zmax);
data(isnan(data)) = 0;
data = smooth3(data,'box',15);

%% 3D-isosurfaces view, Draw
tic
disp('plotting data...'); 



% adust intensity threshold for nascent transcript detection
lev1 = 6;
lev2 = .92; 

figure(4); clf;
if shownucs ==1
patch(isosurface(Y,X,Z,data,nt,'material','dull'),'FaceColor','blue','EdgeColor','none'); hold on;  alpha .3; 
end

if disp_nascent == 1
patch(isosurface(Y,X,Z,Isect1(1:stp:end,1:stp:end,1:Zmax),lev1*25000,'material','shiny'),'FaceColor','red','EdgeColor','none'); alpha .3
patch(isosurface(Y,X,Z,Isect2(1:stp:end,1:stp:end,1:Zmax),lev2*15100,'material','shiny'),'FaceColor','green','EdgeColor','none');  
end

set(gcf,'color','k'); set(gca,'color','k');


toc

%%  3D-isosurfaces Draw dots
tic

points1 = [d1(:,1)*50-yp1*50,d1(:,2)*50-xp1*50,zspace*Zs-zspace*d1(:,3)]; 
points2 = [d2(:,1)*50-yp1*50,d2(:,2)*50-xp1*50,zspace*Zs-zspace*d2(:,3)]; 

V1 = zeros(size(X));  V2 = V1; 

ind1 = cell(Zs,1);
ind2 = cell(Zs,1);
for zs =1:Zs; 
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
V1(inds1) = 20;
V2(inds2) = 20;


figure(4); hold on;
 patch(isosurface(X,Y,Z,V1,1,'material','dull'),'FaceColor','c','EdgeColor','none');
patch(isosurface(X,Y,Z,V2,.5,'material','dull'),'FaceColor','k','EdgeColor','none');
% camlight left; 
%light('Position',[1000,1.5E4 8000],'Style','local');
% lighting phong; % view(112,34); 
view(194,80);
%  

toc


