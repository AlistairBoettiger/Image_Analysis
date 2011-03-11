%%                      Unsupervised_DotFinding.m
%
% Alistair Boettiger                                   Date Begun: 03/10/11
% Levine Lab                                        Last Modified: 03/10/11
%

clear all;


% Input options 
rawfolder = '/Volumes/Data/Lab Data/Raw_Data/02-06-11/';
stackfolder = 'MP10_22C_sna_y_c/';
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
fname ='MP10_22C_sna_y_c';
%fname = 'MP05_22C_sna_y';

getpreciseZ = 0;
consec_layers = 2;

   plotdata = 0 ;% don't show 
   showhist = 1; % show histogram of mRNA counts per cell. 
   bins = 40; 
   mRNA_channels =  2; % 1; % total mRNA channels

    showim = 1; % show colorcoded mRNA counts per cell
    t = .2; % threshold for region definition plotting
    spread = 1.5; % over/under

%---- Dot Finding Parameters ----- %
    sigmaE = 3;%  IMPORTANT
    sigmaI = 4; % IMPORTANT
    min_int  = 0.02;% .05 % not necessary Fix at Zero
    FiltSize = 30;% 
    min_size = 10;% 
   
    % Build the Gaussian Filter   
    Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
    Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
    Filt = Ex -Ix;
%---------------------------------%

Data = cell(1,mRNA_channels); 

for e=2:2 % 100
    tic 
    disp('loading data...');
    if e<10
        emb = ['0',num2str(e)];
    else
        emb = num2str(e);
    end
   
    try load([folder,fname,'_',emb,'_nucdata.mat']); 
            %
            % Loads the following variables.
            %     NucLabeled = downscaled labeled map 
            %     nuc_cents = nuclei centroids
            %     Nucs = downsized raw nuclei image % NOT SAVED
            %     conn_map = connectivity matrix
            %     Cell_bnd = image map of cell boundaries 
    catch me
        disp(me.message)
        break
    end
       
    filename = [rawfolder,'/',fname];
    Im = lsm_read_mod([filename,'.mat'],str2double(emb),1.5E4);    
    
% Focus on subset of image: 
     m =    .9; % 1/2048;  % .7; %   1/2048; %  .85;
    Zs = length(Im);
    [h,w] = size(Im{1,1}{1}); 
    xp1= floor(h/2*m)+1; 
    xp2 = floor(h/2*(2-m))+1;
    yp1 = floor(w/2*m)+1;
    yp2 = floor(w/2*(2-m))+1;
    hs = yp2-yp1+1; 
    ws = xp2-xp1+1;
    disp(['Coordinates:  ', num2str(xp1), ' : ', num2str(xp2), ',   ' num2str(yp1), ' : ', num2str(yp2) ] );
   
    toc
    
    for mRNAchn = 1:mRNA_channels
            DotData = cell(1,Zs);    
            DotMasks = cell(1,Zs); 
            tic; disp('finding dots...'); 
            for z = 1:Zs % z = 20         
                  [cent1,bw1,dL] = dotfinder(Im{1,z}{mRNAchn}( xp1:xp2,yp1:yp2 ),Ex,Ix,min_int,min_size);
                  DotData{z} = cent1;
                  DotMasks{z} = dL;
            end     
            toc;
            
            Cents = cell2mat(DotData');
            
        %%
        % consec_layers = 2;
        
        dotC = CheckDotUpDown(DotData,DotMasks,Im,mRNAchn,hs,ws,plotdata,getpreciseZ,consec_layers);
        % Project all layers
         Imax = imread([rawfolder,stackfolder,fname,'_',emb,'_max.tif']); 
          Imax_dots = Imax(xp1:xp2,yp1:yp2,mRNAchn); 
         %Imax_dots = Imax(:,:,mRNAchn); 
         Iout = figure(2); % clf;  imagesc(Imax_dots);
         colordef black; set(gcf,'color','k'); 
        
         colormap hot; hold on;
           % plot(  xp1+dotC(:,1),yp1+dotC(:,2),'bo','MarkerSize',4 );
            plot(  dotC(:,1),dotC(:,2),'w+','MarkerSize',14 );
            plot(  Cents(:,1),Cents(:,2),'yo','MarkerSize',4);
          %  saveas(Iout,[folder,fname,'_chn',num2str(mRNAchn),'.fig']); 
        %%
         inds = floor(dotC(:,2))+floor(dotC(:,1))*h;   
         inds(inds>w*h) = w*h;        
        
  % % $$$$$$$ % Loop through nuclei counting total dots in region % $$$$$$$$$ % %        
         Nnucs = max(NucLabeled(:)); % total nuclei 
         hn = size(NucLabeled,1);  % size of rescaled nuclear image
         NucLabeled = imresize(NucLabeled,h/hn,'nearest'); % upscale NucLabeled to resolution of mRNA chanel;  
                  
        % Get list of all pixels associated with each nucleus               
        imdata2 = regionprops(NucLabeled,'PixelIdxList','Area'); 
        C=NucLabeled;
        mRNA_cnt = zeros(1,Nnucs); % store counts of mRNA per cell  
        mRNA_den = zeros(1,Nnucs);  % store densities of mRNA per cell
        nuc_area = zeros(1,Nnucs); 
        for i=1:Nnucs
            mRNA_cnt(i) = length(intersect(imdata2(i).PixelIdxList,inds));
            C(C==i) = mRNA_cnt(i);
            mRNA_den(i) = mRNA_cnt(i)/length(imdata2(i).PixelIdxList)  ; 
            nuc_area(i) = length(imdata2(i).PixelIdxList);
        end
        % normalize density to the average cell area
        mRNA_sadj = mRNA_den*mean([imdata2.Area]);
        
        % more stats  
        m_cnt = mean(mRNA_cnt);
        s_cnt = std(mRNA_cnt);
        m_den = mean(mRNA_sadj);
        s_den = std(mRNA_sadj);  
            % save([handles.fdata,'/','test']);
            % load([handles.fdata,'/','test']);    
            
      %% Plotting counts 
         if showhist == 1
                colordef white; 
                figure(5); clf; hist(mRNA_cnt,bins); set(gcf,'color','w');
                title(['mRNA per cell. mean = ',num2str(m_cnt,4),' std=',num2str(s_cnt,4)]); 
                figure(4); clf; hist(mRNA_sadj,bins);set(gcf,'color','w');
                title(['Cell size adjusted mRNA per cell. mean = ',...
                num2str(m_den,4),' std=',num2str(s_den,4)]); 
            % write to disk? 
         end

         if showim == 1        
            figure(3); clf;  colordef black;
            imagesc(C); colormap('hot'); colorbar; 
            set(gcf,'color','k');  
         end
         
        if t ~= 0 && showim == 1 
            Fig_regvar = figure(40); subplot(1,2,mRNAchn);
            [on_cnts,off_cnts]= fxn_regionvar(NucLabeled,C,mRNA_sadj,t,spread,Nnucs);
        end
    
     %% Export data
     Data{mRNAchn}.nucarea = nuc_area;
     Data{mRNAchn}.dotC = dotC;
     Data{mRNAchn}.mRNAcnt = mRNAcnt;
     Data{mRNAchn}.mRNAden = mRNAden;
     Data{mRNAchn}.mRNAsadj = mRNAsadj;
    end % end loop over mNRA channels
end % end loop over embryos 
