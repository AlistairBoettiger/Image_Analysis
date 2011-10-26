%%                      Unsupervised_DotFinding.m
%
% Alistair Boettiger                                   Date Begun: 03/10/11
% Levine Lab                                        Last Modified: 07/07/11
%

clear all;

tot_time = tic;
% Input options 
old_lab = 0;  Es = 0;  ver = ''; %   '_v4'; % '_v3';% '_v2';
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-06-20/'; %   2011-05-22/'; % 2011-04_and_earlier/'; % % 2011-05-22/'; % 2011-06-20/'; %   '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
rawfolder = '/Volumes/Data/Lab Data/Raw_Data/2011-06-20/'; %  2011-05-22/'; %2011-04_and_earlier/'; %'; % 2011-06-20/'; %  '/Volumes/Data/Lab Data/Raw_Data/02-17-11/'; %%   %

stackfolder =  's07_MP05Hz/';% 's04_MP10/';%   'MP07Hz/'; % 's07_MP08/'; % 'MP12Hz/'; %    's05_MP06/';% 's02_MP01/';% 's01_MP09/';%   'sna2.8Hz/' ;%'s06_MP10_sna18/'; %'s21_MP07/';% 'MP07Hz/';% 's11_G4B/' %  's06_MP10_sna18/'; % %'s10_bcd1x/';%  's11_bcd6x/'; %'s14_comp_cntrl/'; % 's12_cntrl_2label/'; %'MP02_22C/'; %'MP01_22C/'; % 'MGa1x/'; % 'MP10_22C/'; %'MP05_22C/'; %'YW_ths_sog/'; % 'MP10_22C/'; %  % 'MP09_22C/'; % 'MGa2x/'; % 'MGa1x/'; % 'MGa2x/'; % 'MP10_22C_sna_y_c/'; %
fname ='s07_MP05Hz_22C'; ver = '_v2'; % 'MP10Hz_c'; %'MP07Hz_snaD_22C_b' ; ver = '_v3';%  'MP08Hz_snaD_22C_b'; % 'MP12Hz_snaD_22C_b'; %    's04_MP10Hz'; % 's05_MP06Hz_b';% 's02_MP01_Hz_22C_b'; % 's01_MP09_Hz_22C_c'; %'sna2.8Hz_snaD_22C'; % 's06_MP10_sna18_b'; % 'MP07het_snaD_22C'; %  'MP07Hz_snaD_22C';%'s11_G4B_LacZ';% 's06_MP10_sna18_b'; % 's05_MP06Hz'; %   %'s10_bcd1x';% 's11_bcd6x'; % 's14_comp_cntrl'; Es =1; % 's12_cntrl_2label'; Es = 1; % 'MP09_22C_hb_y_f'; Es = 7; %  'MP02_22C_hb_y'; Es = 9; % 'MP02_22C_hb_y_b'; Es = 10; %  % 'MP01_22C_hb_y_f'; Es = 12; % 'MP01_22C_hb_y_c'; Es = 10; % 'MP01_22C_hb_y'; Es = 13; % 'MGa1x_LacZ_b'; Es = 12; %  'MP10_22C_sna_y_e'; Es = 12; %  'MP05_22C_sna_y_c'; Es =7; %  'MP10_22C_sna_y_d3'; Es = 1;  %'YW_ths_sog'; Es = 12;  % % 'MP09_22C_hb_y_e'; Es = 10; % 'MP09_22C_hb_y_d'; Es=11; % 'MGa2x_LacZ_sna_b'; Es = 10; % 'MP10_22C_sna_y_d';   % 'MGa_LacZ'; %'MGa2x_LacZ_sna'; %'MP10_22C_sna_y_c'; old_lab = 1;  % 'MP05_22C_sna_y'; old_lab = 1; % 
mRNA_channels = 2;% 2; %  3; %  1; % total mRNA channels

sname =  fname; % 'MP07het_snaD_22C_1';% '_1'; % additional label on slide. 


% MP10_22C_sna_y_c and MP05_22C all done at 3.5, 4, 0.03, 30, 30
% MGa2x and MGa1x all done at 2.5, 3, 0.03, 30, 30



%%

filename = [rawfolder,'/',fname];     
load([rawfolder,stackfolder,sname,'.mat'])  

w = Datas.Stack1.Image1.IMG.width;
h = Datas.Stack1.Image1.IMG.height; 
if Es==0
    Zs = Datas.LSM_info.DimensionZ; 
    Es = length(fields(Datas)) - 3;   % Number of Stacks
end
% ------- Option: Focus on subset of image: ------------------- %
%      m = 1/2048;  % .9;%    .7; % .5; .7; %   1/2048; % 
%    xp1= floor(h/2*m)+1; xp2 = floor(h/2*(2-m))+1;  yp1 = floor(w/2*m)+1;  yp2 = floor(w/2*(2-m))+1;
%    hs = yp2-yp1+1;     ws = xp2-xp1+1;

ws = 2048; hs = 2048; xp1 = 1; yp1 = 1;
xp2 = xp1 + ws -1; yp2 = yp1 + hs - 1; 
disp(['Coordinates:  ', num2str(xp1), ' : ', num2str(xp2), ',   ' num2str(yp1), ' : ', num2str(yp2) ] );
% ------------------------------------------------------------- %     
    

% -------------- Graphing and Display Options ------------------ %
   show_projected = 0; % show max-project with all dots and linked dots. 
   plotdata = 0; % CheckDotUpDown display parameter
   plotZdata = 0 ;% show z-map of data
   showhist = 1; % show histogram of mRNA counts per cell. 
   showim = 1; % show colorcoded mRNA counts per cell
   bins = 40; % bins for histograms of mRNA
   t = 0; %.45; % threshold for region definition plotting
   spread = 1.3; % over/under
% ------------------------------------------------------------- % 


%---- Dot Finding Parameters ----- %
   % dotfinder's parameters 
    sigmaE = 3;%  IMPORTANT
    sigmaI = 4; % IMPORTANT
  %  min_int  = 0.04;    %  5    ;% .05 % not necessary Fix at Zero
    FiltSize = 30;% 
    min_size = 30;% 
   
  % sphere finding parameters
   getpreciseZ = 0;
  % consec_layers = 3;
   ovlap = 4; 
%---------------------------------%



    % Build the Gaussian Filter   
    Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
    Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
 
    disp('Running DotFinder1'); 
  
%%
for e= 1:Es
%%
disp('loading data...');
    tic 
  
    if e<10
        emb = ['0',num2str(e)];
    else
        emb = num2str(e);
    end
    disp(['analyzing embryo, ',emb,'...']);
   
    
    try load([rawfolder,stackfolder,fname,'_',emb,'_nucdata.mat']);
    
    catch err
        disp(err.message)
        try load([folder,fname,'_',emb,'_nucdata.mat']);
            
            % Loads the following variables.
            %     NucLabeled = downscaled labeled map 
            %     nuc_cents = nuclei centroids
            %     Nucs = downsized raw nuclei image % NOT SAVED
            %     conn_map = connectivity matrix
            %     Cell_bnd = image map of cell boundaries 
        catch me
            disp(me.message)
            disp('trying next embryo...'); 
             continue
        end
    end
   
    toc
    
    
    for mRNAchn = 1:mRNA_channels % mRNAchn =2
        
         if mRNAchn == 1;
               min_int  = 0.05;  % just for speed 
               consec_layers = 3;
         elseif mRNAchn == 2
               min_int  = 0.02; %
                consec_layers = 3;
         elseif mRNAchn == 3
               min_int = .01;
               consec_layers = 3;
         end
        
            DotLabels= cell(1,Zs); 
            DotData = cell(1,Zs);    
            Inds = cell(1,Zs); 
            Ints = cell(1,Zs); 
            im_folder = cell(1,Zs);
            
            tic; disp('finding dots...'); 
            for z = 1:Zs % z = 11     
                 im_folder{z} = [rawfolder,stackfolder,fname,'_',emb,'_z',num2str(z),'.tif'];
                 try
                 Iin_z = imreadfast(im_folder{z});       
                 [DotLabels{z},DotData{z},Inds{z},Ints{z}]  = dotfinder(Iin_z(xp1:xp2,yp1:yp2,mRNAchn),Ex,Ix,min_int,min_size);
                 catch err
                     disp(err.message); 
                     Zs = z-1; 
                    break
                 end 
            end
            toc;
            
            % resize;        
            DotLabels= DotLabels(1:Zs); 
            DotData = DotData(1:Zs);    
            Inds = Inds(1:Zs); 
            Ints = Ints(1:Zs); 

            
        %%

         intype = class(Iin_z);
         dotC =  CheckDotUpDown(DotLabels,DotData,Inds,Ints,plotdata,getpreciseZ,consec_layers,ovlap,xp1,xp2,yp1,yp2,intype);
         Cents = cell2mat(DotData');
         
        % Project all layers
         
        if show_projected == 1
            try
                Imax = imread([rawfolder,stackfolder,fname,'_',emb,'_max.tif']); 
            catch err
                disp(err.message); 
                Imax = imread([rawfolder,stackfolder,'max_',fname,'_',emb,'.tif']); 
            end
            
            Imax_dots = Imax(xp1:xp2,yp1:yp2,mRNAchn);  
            figure(2); 
            Iout = figure(2);  clf;  imagesc(Imax_dots);
            colordef black; set(gcf,'color','k'); 
            colormap hot; hold on;
            plot(  dotC(:,1),dotC(:,2),'w+','MarkerSize',14 );
            plot(  Cents(:,1),Cents(:,2),'yo','MarkerSize',4);
            saveas(Iout,[folder,fname,'_',emb,'_chn',num2str(mRNAchn),ver,'.fig']); 
        end
        %%
        
    clear Imax Cents DotData DotLabels Inds Ints Iin_z 
        
        %%
        
        tic
        disp('assigning dots to nuclei...');
        inds = floor(dotC(:,2))+floor(dotC(:,1))*hs;   
        inds(inds>ws*hs) = ws*hs;        
        
  % % $$$$$$$ % Loop through nuclei counting total dots in region % $$$$$$$$$ % %        
        
         hn = size(NucLabeled,1);  % size of rescaled nuclear image
         
         NucLabel = imresize(NucLabeled,h/hn,'nearest'); % upscale NucLabeled to resolution of mRNA chanel;  
         NucLabel = NucLabel(xp1:xp2,yp1:yp2);
         %figure(3); clf; imagesc(NucLabel);
         Nend = max(NucLabel(:)); % total nuclei 
          
          Nmin = single(NucLabel); 
          Nmin(Nmin==0)=NaN; 
          Nstart = min(Nmin(:)); 
          Nucs_list = nonzeros(unique(NucLabel));
          Nnucs = length(Nucs_list);
          
%           M = NucLabel;
%           M(inds) = 300; 
%           figure(1); clf; imagesc(M);
         
   %     % Get list of all pixels associated with each nucleus         
   %     % imdata2 = regionprops(NucLabeled,'PixelIdxList','Area'); 
        
     %   C=NucLabel;
        mRNA_cnt = zeros(1,Nnucs); % store counts of mRNA per cell  
        mRNA_den = zeros(1,Nnucs);  % store densities of mRNA per cell
        nuc_area = zeros(1,Nnucs); 
        if showim == 1
            Plot_mRNA = single(NucLabeled);
        end
        for i=1:Nnucs; % i = 4
            nn = Nucs_list(i);
            imdata.Area(i) = length(find(NucLabel==nn));
            imdata.PixelID{i} = find(NucLabel==nn);
            mRNA_cnt(i) = length(intersect(imdata.PixelID{i},inds));
         %   C(NucLabel==nn) = mRNA_cnt(i);
            mRNA_den(i) = mRNA_cnt(i)/imdata.Area(i); 
            nuc_area(i) = length(imdata.PixelID{i});
            if showim == 1
                Plot_mRNA(NucLabeled==nn) = single(mRNA_den(i));
            end
        end
        % normalize density to the average cell area
        mRNA_sadj = mRNA_den*mean(imdata.Area);
        
        % more stats  
        m_cnt = mean(mRNA_cnt);
        s_cnt = std(mRNA_cnt);
        m_den = mean(mRNA_sadj);
        s_den = std(mRNA_sadj);  
            % save([handles.fdata,'/','test']);
            % load([handles.fdata,'/','test']);    
        toc
            
            
      %% Plotting counts 
      tic
      disp('plotting and saving data...');
         if showhist == 1
                colordef white; 
%                 figure(5); clf; hist(mRNA_cnt,bins); set(gcf,'color','w');
%                 title(['mRNA per cell. mean = ',num2str(m_cnt,4),' std=',num2str(s_cnt,4)]); 
                histfig  = figure(25); clf; 
                hist(mRNA_sadj,bins);
                set(gcf,'color','w');
                title(['Cell size adjusted mRNA per cell. mean = ',...
                num2str(m_den,4),' std=',num2str(s_den,4)]); 
%            saveas(histfig,[folder,fname,'_',emb,'_chn',num2str(mRNAchn),'_hist',ver,'.jpg'],'jpg'); 
            % write to disk? 
         end

         if showim == 1        
            mRNA_map = figure(3); clf;  colordef black;
            imagesc(Plot_mRNA*mean(imdata.Area)); colormap('hot'); colorbar; 
            set(gcf,'color','k');  
            saveas(mRNA_map,[folder,fname,'_',emb,'_chn',num2str(mRNAchn),'rvar',ver,'.jpg'],'jpg'); 
         end
         
         
        if t ~= 0 && showim == 1 
            Fig_regvar = figure(40); clf; % subplot(1,2,mRNAchn);
            [on_cnts,off_cnts]= fxn_regionvar(NucLabel,Plot_mRNA,mRNA_sadj,t,spread,Nnucs,Nucs_list);
             saveas(Fig_regvar,[folder,fname,'_',emb,'_chn',num2str(mRNAchn),'rvar',ver,'.fig']); 
        end
    
     clear imdata M C W  mRNA_map Fig_regvar histfig Iout  
        %
     %% Export data

Rpars.sigmaE = sigmaE;
Rpars.sigmaI = sigmaI;
Rpars.min_int = min_int;
Rpars.FiltSize = FiltSize;
Rpars.min_size = min_size;
Rpars.getpreciseZ = getpreciseZ;
Rpars.consec_layers = consec_layers;
Rpars.ovlap = ovlap; 
     
       save([folder,fname,'_',emb,'_chn',num2str(mRNAchn),'_data',ver,'.mat'],...
           'nuc_area','dotC','mRNA_cnt','Plot_mRNA','mRNA_sadj','Rpars'); 
     
     clear nuc_area dotC mRNA_cnt mRNA_sadj Plot_mRNA 
    
     toc
    end % end loop over mNRA channels
       %  clean up;
        clear Iin_z DotData DotMasks I_max cent1 bw dL Cents ...
            Nmin imdata imdata2 NucLabel NucLabeled Plot_mRNA M C ...
            nuc_area dotC mRNA_cnt mRNA_den mRNA_sadj inds conn_map;
            
    
end % end loop over embryos 


      clear Iin_z DotData DotMasks I_max cent1 bw dL Cents ...
            Nmin imdata imdata2 NucLabeled Plot_mRNA M C ...
            nuc_area dotC mRNA_cnt mRNA_den mRNA_sadj;
              
       Tout = toc(tot_time)/(60*60);
      disp(['elpased time = ',num2str(Tout), ' hours']); 
      
      disp('All slide data saved'); 
      
        

