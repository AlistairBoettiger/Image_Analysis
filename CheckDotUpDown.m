%%                              CheckDotUpDown
% Alistair Boettiger                                   Date Begun: 01/30/11
% Levine Lab                                        Last Modified: 11/22/11


%% Description
% Stitch multi-stack data together to find 3D positions of all dots. 
% avoids multiple dot counting and avoids fusion of dots in Z.

% if getpreciseZ is off, the first layer in which the dot occurs is used as
% the z postion, rather than the brightest layer  

function [New_dotC,Linx,Liny] = CheckDotUpDown(DotLabels,DotData,Inds,Ints,plotdata,getpreciseZ,consec_layers,ovlap,xp1,xp2,yp1,yp2,intype,watershedZ)

%% Updates
% Rewritten 03/07/11 to convert more things to uint16 / uint8 to save
% memory (even fragment of single stack is several gigs of active mem). 
%  Reworked 06/06/11 to use only image data computed in dot finder, no
%  action on raw image data.  Speeds up incredibly.  
%

%% Approach 3: 
% use 'ovlap'-size pixel masks instead of min distance for speed
% use linear indexing of all dots
tic 

%  plotdata = 1;  h = hs; w = ws; getpreciseZ =1; ovlap = 2; intype = 'uint16';  
% DotData = DotData1; Inds = Inds1; Ints = Ints1; DotLabels = DotLabels1; 
% DotData1 = DotData; Inds1 = Inds; Ints1 = Ints; DotLabels1 = DotLabels; 
 % DotData = DotData(mRNAchn,:); Inds = Inds(mRNAchn,:);  Ints = Ints(mRNAchn,:); DotLabels = DotLabels(mRNAchn,:);
 % mRNAchn = 2; 
 
 comp_onVoff = 0;
hs = yp2 - yp1+1; 
ws = xp2 - xp1+1;

% create list of 3d corrdinates of all dots.  Also assigns all dots a
% unique linear index.
Zs = length(DotData);
dotsinlayer = zeros(1,Zs,'single'); % important for later math functions
dotzpos = cell(Zs,1); 
for z = 1:Zs
    dotsinlayer(z) = size(DotData{z},1);
    dotzpos{z} = z*ones(dotsinlayer(z),1);
end
dotC = [cell2mat(DotData'), cell2mat(dotzpos)];  
%clear DotData; 


NDots = length(dotC); % total number of dots;
maxdots = cast(max(dotsinlayer) +100, 'uint16');
disp(['Max dots per layer = ',num2str(maxdots)]); 
disp(['Total dots = ',num2str(NDots)]); 

        % Rather memory inefficient, I have the same centroid data stored
        % in 2 different data structures.  Could build this guy to start
        % with.  

%% Troubleshooting plots (3D dot scatter)       
% figure(10); clf;  plot3(dotC(:,1),dotC(:,2),dotC(:,3),'w.');
% % for n=1:NDots
% %     z = dotC(n,3); 
% %   figure(10); hold on;
% %   plot3(dotC(n,1),dotC(n,2),dotC(n,3),'color',[0,1,0],'Marker','.','MarkerSize',10);
% % end
% %       
%  hold on; 
%   for n = 1:NDots  % n = 300
%       z = dotC(n,3); 
%       clear Iw;  
%       Iw = imread(im_folder{z}); 
%        % figure(4); clf; imagesc(3*Iw(:,:,mRNAchn));
%         Iw = Iw(xp1:xp2,yp1:yp2,mRNAchn);
% %         figure(4); clf; imagesc(Iw); hold on; 
% %         plot(round(dotC(n,1)),round(dotC(n,2)),'w+');
%         
%        cl = double(3*Iw( round(dotC(n,2)),round(dotC(n,1)),:))/255;
%        figure(10); hold on;
%       plot3(dotC(n,1),dotC(n,2),dotC(n,3),'color',cl.*[0,1,1],'Marker','.','MarkerSize',10);
%      % pause(.001); 
%   end
%  
%  figure(11); clf; plot(DotData{20}(:,1),DotData{20}(:,2),'w.');
%  hold on; plot(DotData{21}(:,1),DotData{21}(:,2),'c.');
% hold on; plot(DotData{22}(:,1),DotData{22}(:,2),'b.');
        
%%
disp('connecting dots in Z...') 
DotConn = zeros(2*NDots,Zs,'single'); % empty connectivity matrix for all dots;  as a uint16 restricts this to 65,536 dots per image.  
ConnInt = zeros(2*NDots,Zs,intype); 
LayerJoin = false(2*NDots,Zs); 
% Only enter data in every other line. This leaves black space to allow
% for image segmentation routines to be used that will treat each dot as
% separate.  

% pre-calc
Rs = cell(Zs,1);
for z=1:Zs
         Rz = zeros(hs,ws,'uint16'); % this is uint16 just for data size maintance.    
         Rz(Inds{z}) = DotLabels{z}; % convert indices to raster map 
         Rs{z} = imdilate(Rz,strel('disk',ovlap));
end


for Z = 1:Zs % The primary layer Z = 8
         inds1 = Inds{Z};      
         st1_dot_num = sum(dotsinlayer(1:Z-1)); % starting dot number for the layer under study     
               
    for z=1:Zs % compare primary layer to all other layers  z = Z+1 
      % Need to get linear index to stick correctly in array of all dots.  
         stz_dot_num = sum(dotsinlayer(1:z-1));  % starting dot number for the comparison layer     
         inds_zin1 = Rs{z}(inds1);% dot-indicies of layer Z overlapping layer z.       
         indsT = single(inds_zin1) + stz_dot_num; % convert layer specific dot-indices to total overall dot-indices 
         indsT(indsT == stz_dot_num) = 0; % makes sure missing indices are still 'missing' and not last of previous layer.   
         
%         % trouble shooting        
%          figure(1); clf; 
%             im_folder{Z} = [rawfolder,stackfolder,fname,'_',emb,'_z',num2str(Z),'.tif'];
%             Iin_z = imreadfast(im_folder{Z});  
%             I_test = Iin_z(xp1:xp2,yp1:yp2,mRNAchn+1);  imagesc(I_test); colorbar; colormap hot;
%              hold on; plot(dotC(indsT,1),dotC(indsT,2),'c+');
         
         DotConn(2*st1_dot_num+1:2:2*(st1_dot_num + dotsinlayer(Z)),z) =  indsT; % STORE in DotConn matrix the indices 
         ConnInt(2*st1_dot_num+1:2:2*(st1_dot_num + dotsinlayer(Z)),z) = Ints{z}(inds_zin1+1); % Ints{z}(inds_zin1+1); ; %  % store actual intensities.  
        %  figure(3); clf; imagesc(DotConn); shading flat;
    end
    LayerJoin( 2*st1_dot_num+1 :2*(st1_dot_num + dotsinlayer(Z)),Z) = true(2*dotsinlayer(Z),1);    
end
toc


% figure(3); clf; imagesc(ConnInt); colormap hot; colorbar;





%%  Trouble-shooting: Draw lines between connected dots
% 
% 
% figure(1); clf; 
%   Imax = imread([rawfolder,stackfolder,'max_',fname,'_',emb,'.tif']); 
%             Imax_dots = 3*Imax(xp1:xp2,yp1:yp2,1:3);  
%             Imax_dots(:,:,3) = .1*Imax_dots(:,:,3);
%             imagesc(Imax_dots);   hold on;
% plot(dotC(:,1),dotC(:,2),'c+');
% for n=1:2*NDots   
%     linx = dotC(nonzeros(DotConn(n,:)),1);
%     liny = dotC(nonzeros(DotConn(n,:)),2) ;  
%     plot(linx,liny,'w');
% end
%       
%    

%%

clear R1 Rz Rs LoZ Inds;

%%

tic
disp('clustering disks...');
% figure(3); clf; imagesc(LayerJoin); 
% figure(3); clf; imagesc(DotConn); colormap hot; shading flat;
% figure(3); clf; imagesc(ConnInt); colormap hot; shading flat;  

% The road blocks in this process are the imageprocessing algorithms which
% require a labeled matrix of regionprops: watershed, bwareaopen

stp = 1000;  
Nsects =  ceil(2*NDots/stp);  
masked_inds = single(zeros(2*NDots,Zs)); 

if comp_onVoff == 1
    % Array to record intensities of all true dots. Used to compare intensity
    % of kept dots to removed dots
    masked_ints = single(zeros(2*NDots,Zs)); 
end
    
%cent = [];
Cent = cell(Nsects,1); 

for k=1:Nsects
% For each dot in the system, there is a row i, which contains all the dots
% above and below it.  We want to concentrate only on the dots in z~=i
% which are connected to the dot as z=i.  For this we use the LayerJoin
% object.  
     mask = DotConn(1+(k-1)*stp:min(k*stp,2*NDots),:)>0; % where dots exist 
     ConnInt_T = immultiply(ConnInt(1+(k-1)*stp:min(k*stp,2*NDots),:),cast(mask,intype)); % report intensities where dot connections are found 
     
     % Keep only main diagnol:
     MD = imadd(cast(LayerJoin(1+(k-1)*stp:min(k*stp,2*NDots),:),intype),ConnInt_T); 
     % figure(4); clf; imagesc(MD);
     MD = MD>0;   
     MD = bwareaopen(MD,20); % Remove all dots in z not connected to the dot at z=i  
     % figure(4); clf; imagesc(MD);
     ConnInt_T = immultiply(ConnInt_T,cast(MD,intype));
   
    mask = ConnInt_T>0;

    % Remove small dots
    mask = bwareaopen(mask,consec_layers); % figure(3); clf; imagesc(mask);
 % bwareaopen is very expensive for big images    
    
 if watershedZ == 1
    % Watershed to split dots
    % W = ConnInt_T.*uint16(mask.*longDots); figure(4); clf; imagesc(W);
    W = immultiply(ConnInt_T,cast(mask,intype)); 
    % figure(4); clf; imagesc(W);
    W = watershed(max(W(:)) - W);   % This is the Memory Kill Step; 
     % figure(3); clf; imagesc(W); colormap lines;
    mask(W==0) = 0; 
    % figure(4); clf; imagesc(mask);
 end
 
    if getpreciseZ == 1
        labeled = bwlabel(mask);
        R = regionprops(labeled,ConnInt_T,'WeightedCentroid');
        tcent =  reshape([R.WeightedCentroid],2,length(R))';
        Cent{k} = [tcent(:,1),(k-1)*stp + tcent(:,2)];
        
        if plotdata == 1;
            figure(4); clf; 
            colordef black; set(gcf,'color','k'); 
            imagesc( immultiply(ConnInt_T,cast(mask,intype)) ); 
            colormap hot; shading flat;  colorbar;
            ylabel('mRNA index'); xlabel('z-depth'); 
            hold on; plot(tcent(:,1),tcent(:,2),'co'); 
            title('Cross-Section of all dots'); 
        end
    end
   
    % Final output is just DotConn split up by mask
       masked_inds(1+(k-1)*stp:min(k*stp,2*NDots),:) =  mask.*DotConn(1+(k-1)*stp:min(k*stp,2*NDots),:)  ;
%      
       
       if comp_onVoff == 1 % For comparison of Dot intensities
        masked_ints(1+(k-1)*stp:min(k*stp,2*NDots),:) =  mask.*single(ConnInt(1+(k-1)*stp:min(k*stp,2*NDots),:) );
       end
end

% figure(3); clf; imagesc(masked_inds); colormap hot;

if plotdata == 1
    cent = cell2mat(Cent); 
end

% [L,N] = bwlabel(masked_inds>0);

toc


 %% Troubleshooting
%  figure(2); clf; 
%  imagesc(Imax_dots); colormap hot;  hold on;
%  plot(dotC(:,1),dotC(:,2),'c+');
% for n=1:2:2*NDots  
%     linx = dotC(nonzeros(masked_inds(n,:)),1); % trailing zero prevents cellfun error   
%     liny = dotC(nonzeros(masked_inds(n,:)),2);
%     plot(linx,liny,'c'); 
%      text(dotC((n+1)/2,1)+1,dotC((n+1)/2,2),[' ', num2str(dotC((n+1)/2,3) )],'color','c','FontSize',8);
% end
% 
% 
%  figure(1); clf; 
%  imagesc(Imax_dots(:,:,mRNAchn)); colormap hot;  hold on;
% for n=1:2:2*NDots  
%     linx = dotC(nonzeros(masked_inds(n,:)),1); % trailing zero prevents cellfun error   
%     liny = dotC(nonzeros(masked_inds(n,:)),2);
%     plot(linx,liny,'c'); 
% %     text(dotC((n+1)/2,1)+1,dotC((n+1)/2,2),[' ', num2str(dotC((n+1)/2,3) )],'color','c','FontSize',8);
% end


%%  Free up some Memory
%clear  ConnInt_T DotConn LayerJoin  mask  ConnInt
% 

%% New method for removing duplicates

% masked inds has continuous blocks of connected dots, such that
% dotC(masked_inds) gives the the x and y coordinates of of the N
% associated dots with any given mRNA spot.  

% we first pulll out these separate clusters, then we loop through strings
% of dots, comparing them to all other strings of dots, to remove the ones
% which occur twice. 

tic
disp('Building unique-spheres from cross-section disks...');

Linx = cell(NDots,1); Liny = cell(NDots,1); Lind = cell(NDots,1); 
for n=1:2:2*NDots  
    Linx{(n+1)/2} = [dotC(nonzeros(masked_inds(n,:)),1)',NaN,NaN,NaN]; % trailing zero prevents cellfun error   
    Liny{(n+1)/2} = [dotC(nonzeros(masked_inds(n,:)),2)',NaN,NaN,NaN];
    Lind{(n+1)/2} = (n+1)/2; % [nonzeros(masked_inds(n,:))',0,0,0]; % must be a non-loop way to do this
end
%

 UL = cellfun(@(x,y,z) [x(1:3),y(1:3),z(1)],Linx,Liny,Lind,'UniformOutput',0); % Add index to 
 UL = cell2mat(UL);    % this is infact uniform output, not sure why matlab insists we do it this way 
 
 if isempty(UL) == 0
 
     UL = UL(~isnan(UL(:,1)),:);  % remove all empty values
     UL(isnan(UL(:,3)),3) = 0;

     max_clusters = length(UL); 
     ULcat = [[UL(:,1)+UL(:,4);UL(:,2)+UL(:,5);UL(:,3)+UL(:,6)],[UL(:,7);UL(:,7);UL(:,7)]];
     [~,ui] = unique(ULcat(:,1)); % find unique values
     ui = ui(ui<max_clusters);  % only want indices corresponding to original clusters 

     unique_inds = ULcat(ui,2); % take the cluster indices from the set of unique clusters 
     unique_dotX = Linx(unique_inds);
     unique_dotY = Liny(unique_inds);

    Ndots = length(unique_dotX);
    New_dotC = zeros(Ndots,3); 
    for k =1:Ndots
        phalf = floor(length((unique_dotX{k})-3)/2);
        New_dotC(k,1) = unique_dotX{k}(phalf); % median(unique_dotX{k}); %
       New_dotC(k,2) =   unique_dotY{k}(phalf); %median(unique_dotY{k}); %

    %     [New_dotC(k,1),idx] = mymedian(unique_dotX{k}(1:end-3)); %
    %    % if isempty(idx)~=0
    %         New_dotC(k,2) = unique_dotY{k}(idx(1)); %
    %    % end

        if plotdata == 1 && getpreciseZ == 1
            New_dotC(k,3) = cent(cent(:,2) == 2*unique_inds(k)-1,1);
        else
            New_dotC(k,3) = dotC(unique_inds(k),3) + phalf;
        end
    end

 else
     Ndots = 0;
     New_dotC = [NaN,NaN,NaN]; 
 end

disp([num2str(Ndots),' total spheres found']);

toc
% 
% 
% figure(4); clf; 
% imagesc(Imax_dots);   hold on;
% plot(New_dotC(:,1),New_dotC(:,2),'w.','MarkerSize',30);
% plot(dotC(:,1),dotC(:,2),'m+');
% 
%    for k=1:length(Linx)
%                 plot(Linx{k}(1:end-3),Liny{k}(1:end-3),'c');
%    end











%% Old Method to Remove duplicates
% 
% tic
% disp('Building spheres from cross-section disks...');
% 
% % figure(3); clf; imagesc(masked_inds); colormap hot;
% 
% % testL = NaN*zeros(1,2*NDots); 
% remove_dot = zeros(NDots,1); 
% stacked_dots =0;
% % loop through all dots
% 
% for i = 1:2:2*NDots-1 % i =15 
%     j = find(masked_inds(i,:)); % non-zero indices of masked_inds matrix at i
%     counted = masked_inds(i,j(2:end));  % global-dot-indices of overlapping dots. 
%     if isempty(j) == 0
%         stacked_dots = max(j)-min(j) > length(j)-1;
%     
%         if stacked_dots == 0  && getpreciseZ == 1
%              ii = find(cent(:,2)==i,1);
%              dotC((i+1)/2,3) = cent(ii(1),1);
%              i_self = masked_inds(i,j(1));
%              i_max = (cent(ii(1),2)+1)/2;
%              dotC(i_self,1) = dotC(i_max,1); 
%              dotC(i_self,2) = dotC(i_max,2); 
%         end
%     else
%         remove_dot((i+1)/2) = 1;
%     end
%     if stacked_dots == 1% if stacked dots split up.  
%         brk_pts =[0, find(diff(j)>1),length(j),length(j)]; % breakpoints in stack 
%         possibles = masked_inds(i,j); % all possible multicounted indices 
%         ii = find(possibles == (i+1)/2,1) ; % find this breakpoint    
%         % only need this if low intensity points have been removed
%         if isempty(ii); [jnk, ii] = min( ((i+1)/2 - possibles).^2 );  end
%        % find nearest breakpoint without going over
%           kk = (ii-brk_pts); kk(kk<0) = 100; [jnk,bi] = min(kk);    
%           counted = possibles(brk_pts(bi)+2:brk_pts(bi+1)); 
%           if  getpreciseZ == 1
%               ii = find(cent(:,2)==i);
%               ii = ii( min(bi,length(ii)) );
%               dotC((i+1)/2,3) = cent(ii,1);
%               i_self = masked_inds(i,j(1));
%               i_max = (cent(ii,2)+1)/2;
%               dotC(i_self,1) = dotC(i_max,1); 
%               dotC(i_self,2) = dotC(i_max,2);         
%           end
%     stacked_dots =0;      
%     end
%      remove_dot(counted) = 1;   
%    % testL(i) = length(remove_dot);
% end
% toc
% sum(remove_dot);
% 
% % clear DotData DotMasks 
% % clear ConnInt ConnInt_T DotConn  LayerJoin R1 Rz Rs  
% % clear W Nuclabeld  Cell_bnd dL bw1 
% % clear dotC New_dotC test_dotC  Cents
% % clear Iin_z LoZ MD conn_map mask masked_inds 
% 
% % figure(1); clf; plot(testL,'w.');
% % sum(stacked_dots)
%           % NB can't sum stacked dots, all stack dots are also multiply
%           % counted.  i.e. the first time a doublet is enountered we say
%           % Tstacked = Tstacked + 1, and then we enounter that the other of
%           % the pair and again say Tstacked = Tstacked + 1;  
% 
%  if plotdata ==1
%     test_dotC = dotC; test_dotC(logical(remove_dot),3) = 0;
%     figure(5); clf; imagesc(masked_inds); colormap hot;
%     hold on; plot(test_dotC(:,3),linspace(1,2*NDots-1,NDots),'co');
%  end
%  
%  New_dotC = dotC(~remove_dot,:);
% figure(2); clf; 
% imagesc(Imax_dots);   hold on;
% plot(New_dotC(:,1),New_dotC(:,2),'c+');
 
    %%      
          
% %New_dotC = dotC(~remove_dot,:);
% dotC = dotC(~remove_dot,:);
% [N_dots, jnk] = size(dotC); 
% %N_dots = NDots - sum(remove_dot) % sum(stacked_dots)
% %N_dots = length(New_dotC); 
% disp(['Counted ',num2str(N_dots),' spheres']); 
% 
% 
% 
% if comp_onVoff == 1
%     figure(7); clf;  subplot(2,1,1);
%     hist(log(nonzeros(masked_ints(~remove_dot,:)))); title('intensities of kept dots');
%     subplot(2,1,2);
%     hist(log(nonzeros(masked_ints(logical(remove_dot),:)))); title('intensities of removed dots');
% end



% figure(2); clf; 
% imagesc(Imax_dots);   hold on;
% plot(dotC(:,1),dotC(:,2),'c+');



