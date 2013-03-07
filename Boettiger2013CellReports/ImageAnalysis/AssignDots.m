
% Assign dots to nuclei 

function  [mRNA_sadj,Plot_mRNA] = AssignDots(dotC,NucLabeled,imaxes)

showim = 1;
showhist = 1; 
 bins = 40; % for histogram

xp1 = imaxes.xp1;
xp2 = imaxes.xp2;
yp1 = imaxes.yp1;
yp2 = imaxes.yp2;
h = imaxes.h;
w = imaxes.w;

ws = xp2 - xp1 +1;
hs = yp2 - yp1 +1;

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
            Plot_mRNA = single(NucLabel);
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
                Plot_mRNA(NucLabel==nn) = single(mRNA_den(i));
            end
        end
        % normalize density to the average cell area
        mRNA_sadj = mRNA_den*mean(imdata.Area);
        % figure(3); clf; imagesc(Plot_mRNA); colorbar;
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
           
            % write to disk? 
         end

         Plot_mRNA = Plot_mRNA*mean(imdata.Area);
         if showim == 1        
            mRNA_map = figure(3); clf;  colordef black;
            imagesc(Plot_mRNA); colormap('hot'); colorbar; 
            set(gcf,'color','k');  
         end
         
         