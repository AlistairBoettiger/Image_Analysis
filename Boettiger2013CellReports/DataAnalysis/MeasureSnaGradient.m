
%% anlz_sna_gradient_data.m

% Alistair Boettiger                                Date Begun: 02/02/11
% Levine Lab                                    Last Modified: 06/23/11


%% Modifications
% Modified 06/12/11 to scale green channel to correct for missed detection
% rate
% Modified 06/14/11 to save oriented mRNA plots with corrected counts.  
% 6/23/11 converted from hb to sna. 
% 5/28/12 Moved to a copy of script into snail project Code

  %% New data-method
  clear all;
  
  maternal = 0;  ver = '';  Es = 15; cor = 0;   nametype = 1;  st_channel = 0;  legon  = 0; vout ='';  bkg_sub = 1; % defaults
  cbar = 1;  chns_flipped = 0; nuc_chn =3; plotnoise = 0;    chns = 1; min_cnt = 0; % 1 % 2; %
  folder = 'C:\Users/Alistair/My Documents/Projects/mRNA_counting/Data/';
   rawfolder = 'D:/Data/';
   % rawfolder = 'K:\Confocal_fall_2011\' %'G:\Raw_Data\';%
  slidedate =  '2013-01-06_no-dist\'; %  '2012-04-01_MP08/'; % '2011-12\';% '2011-05-22/'; % '2011-06-20/'; %'2012-09-06_MLslides\';% '2012-07-27_MLslides/';%  '2012-08-26_mp07/'; % '2012-09-02_MLslides/';% '2012-08-15_MLslides/';%  '2012-08-26_mp07/'; %'2012-08-20_ML124/'; %  '2011-06-20/';% '2011-05-22/'; %  '2012-08-02_MLslides/';%  '2012-08-20_ML124\'; '2012-07-27_MLslides/'; %'2012-08-15_MLslides/'; %   '2012-03-24_4xsna\';%    '2011-02-17/';% '2011-02-17/';% '2011-11/'; %  '2011-04_and_earlier/';  % '2011-06-20/';%              %  
  subfolder =  ''; %'s07_MP08/' ;%'s07_MP05Hz/'; %'s05_MP06/';%  's04_MP10/'; % 's05_MP06/'; %  's05_MP06/';% 's08_MP06Hz\';%  's04_MP10/'; % 's07_MP05Hz/';% 's21_MP07/'; % 'MP10_22C\'; %  'MP05_22C\';% 's08_MP06_cflip\'; %        'MP12Hz/'; %  's07_MP05Hz/'; % 's04_MP10/'; %  'sna2p8het/'; %   'MP07Hz/';% 's04_MP10/';%   ; %    's07_MP08/'; %  's06_MP10_sna18/'; %    'sna2.8Hz/';%   'MP07Hz/';% 's21_MP07/';  % '' %    's05_MP06/'; % 's11_G4B/'; %  % s04_MP10/';%    
      fname =  'MP08Hz_snaD'; ver = '_v2'; % 's6_MP08_b'; % 's05_MP08_b';        ver = '';  %'MP08Hz_snaD_22C_b'; %  's07_MP05Hz_22C_c'; ver ='_vN3'; %'s05_MP06Hz_b';ver = '_vN4'; chns=2; %'s04_MP10Hz_b'; ver = '_vN2'; chns =2;%  'ml126';%  'ml128counts'; % 'm197a'; ver = '_v2'; %  'm197b'; %   'm109a'; ver = '_v2';% 'ml118'; chns = 2; ver = '_v2';%  'm116b'; ver='_v2';%  'mp07a2'; %  'mp07b'; Es=20; ver = '_v2'; %'ml124B'; Es=9; %  'ml125'; ver = '_v2'; % 'mp10a'; chns = 2; ver = '_vN';%  'm198c'; % % % 's07_MP05Hz_22C'; ver = '_vN2'; chns=2;% 's05_MP06Hz'; ver = '_vN3'; chns =2; %  's08_MP06Hz_cflip'; ver = ''; chns =2; %  'm106a'; Es=16;%   's04_MP10Hz'; ver = '_vN'; chns=2; % 'MP10Hz_c'; '_vN1'; % 's07_MP05Hz_22C'; ver = '_vN1';% 'b308first 2pos'; % 'm198d';%  'm106a';% 'm116b';'m105c'; %  'm105b'; % 'm105a'; % 'm199b_fewissues'; %   'm198d'; %% '4xsna'; %  '4xsna_d'; % 'MP07het_snaD_22C'; %  'snaD_b';  Es = 16; % 'snaD'; % % 's05_MP08_b'; ver = ''; % 's04_MP08'; ver =''; %  's142_sna'; ver = '_v2'; Es = 16; % 'snaD'; % 's140_sna'; ver ='_v3'; % 'MP10_22C_sna_y_e'; ver = '_vN'; % 'MP05_22C_sna_y'; ver = '_vN'; % 'MP10_22C_sna_y_d'; ver = '_vN'; % % 'wt_sna'; ver = '_v2';%   's05_MP06Hz'; ver = '_vN';% 'MP08_snaD_LacZ647'; ver = '_v3';% 'MP05'; ver = '_v3'%   %  %  'MP08Hz_snaD_22C_b' ; %  'MP12Hz_snaD_22C'; ver = '';%'s07_MP05Hz_22C_c';  ver = '';% 'MP10Hz_c'; ver = '_v2'; % 'sna2p8het_30C'; % 'MP07Hz_snaD_22C_b'; ver= '_v3'; %  % 's04_MP10Hz_b'; ver = '_v3'; %  %    's07_MP05Hz_22C_c';%  ver = '_v2';%  'sna2p8het_sim_sna'; ver = '_v2'; Es=1; % 'MP12Hz_snaD_22C_b'; ver = ''; % 's05_MP06Hz_b'; ver = '_v4';%   's06_MP10_sna18_b'; ver = '_v4'; cbar =1;  %     % 's07_MP08Hz_snaD_22C';  % 'sna2.8Hz_snaD_22C';st_channel = 1; % 'MP07Hz_snaD_22C'; % 'MP10_22C_sna_y_d'; ver = '_v3'; % 'MP05_22C_sna_y_c';  ver = ''; %  % 's05_MP06Hz'; ver = '_v2';%   %      's06_MP10_sna18_b'; st_channel = 1;   ver = '_v4'; % ' 's05_MP06Hz_b'; ver = '_v2';  % 's11_G4B_LacZ'; ver = '_v2'; legon =0; %     'MP07het_snaD_22C';%   %  '_v2';  %
  % 
    vc1 = ver; vc2 = ver;   % vc1 = ''; vc2 = '_v4';
  %%
% chns_flipped = 1; 
  
addpath('C:\Users\Alistair\Documents\Projects\GenCode');
addpath('C:\Users\Alistair\Documents\Projects\Image_Analysis');
   vout ='';
  missG = 1.0; 
  ipars{3} = missG; 
  
  manual_orient =[];%  [-10,-45,-35,45,135,45,-145,-50,-80,65];


 
%    try
%      load([folder,fname,'_slidedata',ver], 'Data'); 
%    catch er
%        nametype = 1; 
%    end
  
  
  figure(10); clf; figure(11); clf; figure(12); clf; figure(14); clf;
  
  data = cell(Es,1); 
for e =  1: Es %  e = 5
    if e<10
        emb = ['0',num2str(e)];
    else
        emb = num2str(e);
    end
  
  i = str2double(emb);


      try 
           load([rawfolder,slidedate,subfolder,fname,'_',emb,'_nucdata.mat']);  
          
      catch er
          disp([er.message, ' trying alternate location...']);
           try 
               load([folder,slidedate,fname,'_',emb,'_nucdata.mat']);  
               disp(['found file at ', folder,slidedate,fname,'_',emb,'_nucdata.mat']);
           catch er1
               disp(er1.message);
               disp(['skipping embryo ', emb]);
               continue
           end
      end

      
      if nametype == 1 
          try
               ver = vc1; 
              load([folder,slidedate,fname,'_',emb,'_chn', num2str( st_channel+1),'_data',ver,'.mat']); 
              if chns_flipped == 1  
                    load([folder,slidedate,fname,'_',emb,'_chn2','_data',ver,'.mat']); 
              end
              mRNAsadj = mRNA_sadj; % mRNA_cnt./nuc_area;
             % disp(['chn1 thresh: ', num2str(Rpars.min_int)]);
              ipars{1} = Rpars; 
              if chns ==2
                ver = vc2;
                load([folder,slidedate,fname,'_',emb,'_chn',num2str( st_channel+2),'_data',ver,'.mat']);
                if chns_flipped == 1
                    load([folder,slidedate,fname,'_',emb,'_chn1','_data',ver,'.mat']);
                end
               %  disp(['chn2 thresh: ', num2str(Rpars.min_int)]);
                mRNAsadj2 = mRNA_sadj; % mRNA_cnt./nuc_area;  % 
                ipars{2} = Rpars; 
              end
          catch er
                disp(er.message);
          end
          
      else
        mRNAsadj = Data{i,1}.mRNAsadj;
        mRNAsadj2= Data{i,2}.mRNAsadj;
      end
          
    
%%              
       PlotmRNA = imresize(NucLabeled,.5,'nearest');
       
       if chns == 2; 
            PlotmRNA2 = PlotmRNA; 
       end
       NucLabeled = imresize(NucLabeled,.5,'nearest'); 
                      Nnucs =    max( NucLabeled(:) );
                      for n=1:Nnucs;
                          try
                            PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                          catch er
                              disp(er.message)
                              continue
                          end
                          if chns == 2
                              try
                            PlotmRNA2(PlotmRNA2==n) = mRNAsadj2(n+cor);
                              catch er
                                  disp(er.message);
                              end
                                  
                          end
                      end
   % figure(1); clf; imagesc(PlotmRNA); colormap hot; colorbar;
        
%% Orient image along major axis     

    bw = imresize(makeuint(PlotmRNA,16),.5,'nearest');  
    thresh = graythresh(bw);
    %  figure(3); clf; imagesc(bw);  
    bw = im2bw(bw,thresh);    
   %  figure(3); clf; imagesc(bw);
    bw = imfill(bw,'holes'); 
    bw2 = bwareaopen(bw,1500);
    if sum(bw2(:)) ~= 0;
        bw = bw2;
    end
   % figure(3); clf; imagesc(bw);
  
 
    rprops = regionprops(bwlabel(bw),'MajorAxis','Orientation');
    
 meanvar = 100;
  rotes = 0; 
  %%
 try
 %% 
 while meanvar > 75  && rotes<4;
    
     % automatic rotation by aranging major axis of dominant thresholded object along x-axis.  
     % If this orientation has the max on the left or right extreme variation should be small along the  
   
     
     if isempty(manual_orient)
    NucLabel = imrotate(NucLabeled,(0+90*rotes)-rprops(1).Orientation,'nearest');  
    % figure(1); clf; imagesc(NucLabel);
    rotes = rotes + 1; 
     else
         meanvar = 0;
         NucLabel = imrotate(NucLabeled,manual_orient(e),'nearest'); 
          PlotmRNA_r = imrotate(PlotmRNA,manual_orient(e),'nearest'); 
         % figure(2); clf;  imagesc(PlotmRNA_r); pause(1); 
          if chns == 2;
           PlotmRNA2_r = imrotate(missG*PlotmRNA2,manual_orient(e),'nearest'); 
          end
     end

    
% convert nucleus centroids to indexed postions


    S = regionprops(NucLabel,'Centroid');
    nuc_cents = reshape([S.Centroid],2,length(S));
    [h,w] = size(NucLabel); 
     c_inds = sub2ind([h,w],floor(nuc_cents(2,:)),floor(nuc_cents(1,:)));
    % compute distances
    d = nuc_cents(2,:)*h; % /hn;

    if length(c_inds) < length(mRNAsadj);
        mRNAsadj = mRNAsadj(2:end);
        if chns == 2;
            mRNAsadj2 = mRNAsadj2(2:end);
        end
    end


   % maternal =0; % min(mRNAsadj); %0; % 

     
    %  % Plotting for troubleshooting
    %     C = false(h,w);
    %     C(c_inds) = 1; 
    %     figure(1); clf; imshow(C);

    % Sort by distance from upper left corner (max bcd).  
    c_inds(isnan(c_inds)) = [];
    nuc_order = NucLabel(c_inds);
    [b,m,n] = unique(nuc_order);
    dists = d(m);
    if length(dists) < length(mRNAsadj)
        disp('length dists does not match');
        % dists = [0,dists,dists(end)]; 
       dists = d;
    end
      %  figure(1); clf; plot(dists,mRNAsadj ,'g.');  % check results

    if chns == 2
        Data_notsort = [dists; mRNAsadj ; missG*(mRNAsadj2)]';
    else
        Data_notsort = [dists; mRNAsadj]';
    end
    Data_sort = sortrows(Data_notsort); 
    if isempty(manual_orient)
        meanvar = std(Data_sort(10:20,2));     
         % show rotated image
             PlotmRNA_r = imrotate(PlotmRNA,(0+90*rotes)-rprops(1).Orientation,'nearest');
%          figure(2); clf;  imagesc(PlotmRNA_r)            
        if chns == 2
            PlotmRNA2_r = imrotate(missG*PlotmRNA2,(0+90*rotes)-rprops(1).Orientation ,'nearest'); 
%             figure(3); clf; imagesc(PlotmRNA2_r); colormap hot;           
        end     
    end
 end
 
 % arrange A-P orientation
 PlotmRNA_r = imrotate(PlotmRNA_r,90);
 if chns == 2
     PlotmRNA2_r = imrotate(PlotmRNA2_r,90,'nearest');
 end

 
 
    %%
  catch err
      disp(err.message);
      disp(['error on line ', err.stack.line]);
     %  break
       continue 
 end
 

%% Cluster cells based on distance from anterior pole
Sects = round(sqrt(Nnucs));
Q = cell(1,Sects); 
if chns == 2;
    mu = zeros(Sects,2);
    sigma = zeros(Sects,2);
    bsmu = zeros(Sects,2);
    bssigma = zeros(Sects,2);
    x = zeros(Sects,1); 
    Q2 = cell(1,Sects); 
else
    mu = zeros(Sects,1);
    sigma = zeros(Sects,1);
    bsmu = zeros(Sects,1);
    bssigma = zeros(Sects,1);
    x = zeros(Sects,1); 
end

dbnd = zeros(Sects,1);
for j=1:Sects
    Q{j} = Data_sort( floor((j-1)*Nnucs/Sects) + 1: min(floor(j*Nnucs/Sects),length(Data_sort)), 2);
    
    if chns == 2 
    Q2{j} = Data_sort( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 3);
    mu(j,:) = [nanmean(Q{j}),nanmean(Q2{j})]   ;
    bsmu(j,:) =  [bserr(Q{j},'nanmean'),bserr(Q2{j},'nanmean')] ;
    sigma(j,:) = [nanstd(Q{j}),nanstd(Q2{j})]   ;
    bssigma(j,:) =  [bserr(Q{j},'nanstd'),bserr(Q2{j},'nanstd')] ;
    else
    mu(j,:) = [nanmean(Q{j})]   ;
    bsmu(j,:) =  [bserr(Q{j},'nanmean')] ;
    sigma(j,:) = [nanstd(Q{j})]   ;
    bssigma(j,:) =  [bserr(Q{j},'nanstd')] ;
    end

    x(j) = mean(Data_sort( floor((j-1)*Nnucs/Sects) + 1: min(floor(j*Nnucs/Sects),length(Data_sort)), 1));
    fano = sigma(j)^2/mu(j); 
end

% x = linspace(mean(dists(1:Sects)),mean(dists(end-Sects:1:end)),Sects)*50/1000;

% fix orientation
or = mu(1,1) -  mu(end,1);
if or<0
   Data_sort(:,2) = flipud( Data_sort(:,2));
   mu(:,1) = flipud(mu(:,1));
   sigma(:,1) = flipud(sigma(:,1));
   bssigma(:,1) = flipud(bssigma(:,1)); 
   PlotmRNA_r = fliplr(PlotmRNA_r);

   if chns == 2
        Data_sort(:,3) = flipud( Data_sort(:,3));
        mu(:,2) = flipud(mu(:,2));
        sigma(:,2) = flipud(sigma(:,2));
        bssigma(:,2) = flipud(bssigma(:,2));
        PlotmRNA2_r = fliplr(PlotmRNA2_r); 
   end   
end



%% Plotting


figure(10); subplot(3,ceil(Es/3),e);  colordef white; set(gcf,'color','w');
    plot(Data_sort(:,1),Data_sort(:,2)-bkg_sub*min(mu(:,1)),'k.'); % check results  
    hold on; 
    errorbar(x,mu(:,1)-bkg_sub*min(mu(:,1)),sigma(:,1),'linestyle','none','linewidth',3,'color','r');
    if chns == 2;
        plot(Data_sort(:,1),Data_sort(:,3)-bkg_sub*min(mu(:,2)),'b.'); 
        hold on; 
        errorbar(x,mu(:,2)-bkg_sub*min(mu(:,2)),sigma(:,2),'linestyle','none','linewidth',3,'color','c');
    end
    ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');
    title(['Nuclei = ',num2str(Nnucs)]);
    xlim([min(x),max(x)]); 

    if plotnoise == 1
        figure(11); subplot(3,ceil(Es/3),e);  
        colordef white; set(gcf,'color','w');
        errorbar(x,sigma(:,1)./mu(:,1),bssigma(:,1)./mu(:,1),  'r.','MarkerSize',10); 
        hold on;
        if chns == 2;
            errorbar(x,sigma(:,2)./mu(:,2),bssigma(:,2)./mu(:,2),  'b.','MarkerSize',10); ylim([0,1]);
        end
        plot(x,sqrt(mu(:,1))./mu(:,1),'k.','MarkerSize',10);
        ylabel('CoV'); xlabel('distance (nm)');

        if legon == 1
            if chns == 2
                legend('sna CoV','y CoV','Poisson CoV');
            else
              legend('sna CoV','Poisson CoV');  
            end
        end
        title(['Nuclei = ',num2str(Nnucs)]);
    end

figure(14); subplot(3,ceil(Es/3),e);
 Inuc = imread([rawfolder,slidedate,subfolder,'max_',fname,'_',emb,'.tif']);
 Inuc = Inuc(:,:,nuc_chn);  
 Inuc = imrotate(Inuc,(90+90*rotes)-rprops(1).Orientation,'nearest');
 imagesc(Inuc); colormap hot; %caxis([mean(Inuc(:))/2,2*mean(Inuc(:))]);


 
figure(12); colordef black; subplot(3,ceil(Es/3),e);
imagesc(PlotmRNA_r-bkg_sub*min(mu(:,1))); colormap hot; colordef black;
if cbar == 1
    caxis([min_cnt,max(Data_sort(:,2))-bkg_sub*min(mu(:,1))]);  colorbar;
    
   %  caxis([20,200]); colorbar; axis off;
end
    
if chns == 2
    figure(13); subplot(3,ceil(Es/3),e); 
    set(gcf,'color','k'); colordef black;
    imagesc(PlotmRNA2_r); colormap hot; 
    if cbar == 1 
       % caxis([20,max(Data_sort(:,3))]); colorbar; axis off;
        caxis([min_cnt,200]); colorbar; axis off;
    end
end

    
data{e}.Data_sort = Data_sort;
data{e}.mu = mu;
data{e}.bsmu = bsmu;
data{e}.sigma = sigma;
data{e}.bssigma = bssigma;
data{e}.x = x; 
data{e}.Nnucs = Nnucs; 
data{e}.missG = missG; 
data{e}.ipars = ipars;
data{e}.Inuc = Inuc; 

data{e}.PlotmRNA = PlotmRNA_r; 
if chns == 2
    data{e}.PlotmRNA2 = PlotmRNA2_r;
end
%%
end

figure(12); set(gcf,'color','k');


save([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
% % % 
%   nclass = {'cc14l','cc14l','cc13telo','cc14l'}
% save([folder,slidedate,fname,ver,'_agedata',vout,'.mat'],'nclass'); 
%    disp(['wrote ' folder,slidedate,fname,ver,'_agedata',vout,'.mat']);

%save([folder,slidedate,fname,'_graddata',ver,'.mat'],'data'); 

%%

% e =1;
% figure(1); clf; imagesc(data{e}.PlotmRNA); colordef black; set(gcf,'color','k'); colormap hot; colorbar; caxis([15,110]); axis off;
% figure(2); clf; imagesc(data{2}.PlotmRNA); colordef black; set(gcf,'color','k'); colormap hot; colorbar; caxis([15,300]); axis off;



