%%                              CheckDotUpDown
% Alistair Boettiger                                   Date Begun: 01/30/11
% Levine Lab                                        Last Modified: 02/06/11


%% Description
% Stitch multi-stack data together to find 3D positions of all dots. 
% avoids multiple dot counting and avoids fusion of dots in Z.


function dotC = CheckDotUpDown(DotData,DotMasks,im_folder,mRNAchn,plotdata,getpreciseZ,consec_layers,ovlap,xp1,xp2,yp1,yp2)
%% Updates
% Rewritten 03/07/11 to convert more things to uint16 / uint8 to save
% memory (even fragment of single stack is several gigs of active mem). 
%  
%

%% Approach 3: 
% use 'ovlap'-size pixel masks instead of min distance for speed
% use linear indexing of all dots
tic 

     
intype = 'uint8';

disp('connecting dots in Z...') 
%  plotdata = 1;  h = hs; w = ws; mRNAchn = 1; ovlap = 2;
% DotData = DotData2; DotMasks  = DotMasks2;
  % I = Im{1,z}{mRNAchn}( xp1:xp2,yp1:yp2 );
  
 comp_onVoff = 0;
hs = yp2 - yp1+1; 
ws = xp2 - xp1+1;

% create list of 3d corrdinates of all dots.  Also assigns all dots a
% unique linear index.
Zs = length(DotData);
dotsinlayer = zeros(1,Zs,'single'); % important for later math functions
dotzpos = cell(Zs,1); 
%dotC = [];
for z = 1:Zs
    dotsinlayer(z) = size(DotData{z},1);
    dotzpos{z} = z*ones(dotsinlayer(z),1);
   % dotC = [dotC; DotData{z}, z*ones(dotsinlayer(z),1)];
end
dotC = [cell2mat(DotData'), cell2mat(dotzpos)];  



NDots = length(dotC); % total number of dots;
maxdots = eval([intype,'(max(dotsinlayer) +100)']);
totaldots = length(dotC);
disp(['Max dots per layer = ',num2str(maxdots)]); 
disp(['Total dots = ',num2str(totaldots)]); 

        % Rather memory inefficient, I have the same centroid data stored
        % in 2 different data structures.  Could build this guy to start
        % with.  

%%        
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
%%  
%  figure(11); clf; plot(DotData{20}(:,1),DotData{20}(:,2),'w.');
%  hold on; plot(DotData{21}(:,1),DotData{21}(:,2),'c.');
% hold on; plot(DotData{22}(:,1),DotData{22}(:,2),'b.');
        
   

DotConn = zeros(2*NDots,Zs,'single'); % empty connectivity matrix for all dots;  as a uint16 restricts this to 65,536 dots per image.  
ConnInt = zeros(2*NDots,Zs,intype); 
LayerJoin = false(2*NDots,Zs); 
% Only enter data in every other line. This leaves black space to allow
% for image segmentation routines to be used that will treat each dot as
% separate.  

% pre-calc
Rs = cell(Zs,1);
Inds = cell(Zs,1); 
Ints = cell(Zs,1); 
for z=1:Zs
         inds = floor(DotData{z}(:,2))+floor(DotData{z}(:,1))*hs;  % indices in this layer  
         inds(inds>ws*hs) = ws*hs;  
         Rz = zeros(hs,ws,intype);   
         Rz(inds) = DotMasks{z}(inds); % convert indices to raster map 
         Inds{z} = inds; % pixel-indices of all dots in layer
         Rs{z} = imdilate(Rz,strel('disk',ovlap));
end


% In each layer, index all the intensities at all x-y positions which
% contain a potential dot in any of the z-layers.  This will be a faster
% lookup table for getting the intensity of the pixels at correct
% positions.  
All_pix_inds = cell2mat(Inds);
for z=1:Zs
          clear Iw; 
         Iw = imread(im_folder{z});  % read in whole-image
        Iw = Iw(xp1:xp2,yp1:yp2,mRNAchn); % shrink to working size
        Ints{z} =[0; Iw(All_pix_inds)];
        
end
toc


tic
for Z = 1:Zs % The primary layer Z = 26
         % convert to pixel linear-index
         inds1 = floor(DotData{Z}(:,2))+floor(DotData{Z}(:,1))*hs;  % pixel-indices in this layer  
         inds1(inds1>ws*hs) = ws*hs;  
         R1 = zeros(hs,ws,intype);   
         R1(inds1) = maxdots; % convert indices to raster map   
        % figure(3); clf; imagesc(R1);  colorbar; colormap jet;
        
         st1_dot_num = sum(dotsinlayer(1:Z-1)); % starting dot number for the layer under study     
      
        % figure(3); clf; imagesc(Iw);
         
    for z=1:Zs % compare primary layer to all other layers  z = Z+1 
        clear Loz;
         Loz = R1 + Rs{z}; 
         % figure(3); clf; imagesc(Rs{z});  colorbar; colormap jet;
       %  Loz = R1 + DotMasks{z}; (only makes a small difference to use overlap box rather than within found dot.  
         
%           figure(3); clf; imagesc(Loz);  colorbar; colormap jet;
%          hold on; plot(DotData{Z}(:,1),DotData{Z}(:,2),'w.');
%          hold on; plot(DotData{z}(:,1),DotData{z}(:,2),'c.');
        
        % figure(3); clf; imagesc(DotMasks{z}); 
         Loz(Loz<maxdots+1) = 0; % remove non-overlapping dots;
         Loz = Loz - maxdots; % LoZ is already positive, so we don't need to worry about negative values.  
        %  figure(2); clf;  imagesc(Loz); colorbar; colormap jet;

      % Need to get linear index to stick correctly in array of all dots.  
         stz_dot_num = sum(dotsinlayer(1:z-1));  % starting dot number for the comparison layer     
         inds_zin1 = Loz(inds1); % dot-indicies of layer Z overlapping layer z.       
         indsT = single(inds_zin1) + stz_dot_num; % convert layer specific dot-indices to total overall dot-indices 
         indsT(indsT == stz_dot_num) = 0; % makes sure missing indices are still 'missing' and not last of previous layer.   
         indsI = indsT; indsI(indsI==0) = 1;
         
%          % trouble shooting
%          figure(3); hold on; indsT_NaN = nonzeros(indsT); 
%          for kk = 1:length(indsT_NaN)
%             plot(dotC(indsT_NaN(kk),1),dotC(indsT_NaN(kk),2),'r+');
%          end
         
         DotConn(2*st1_dot_num+1:2:2*(st1_dot_num + dotsinlayer(Z)),z) =  indsT; % STORE in DotConn matrix the indices 
         ConnInt(2*st1_dot_num+1:2:2*(st1_dot_num + dotsinlayer(Z)),z) = Ints{z}(indsI); % Iw(inds_zin1+1) ; %  % store actual intensities.  
        % figure(3); clf; imagesc(DotConn); shading flat;
    end
    LayerJoin( 2*st1_dot_num+1 :2*(st1_dot_num + dotsinlayer(Z)),Z) = true(2*dotsinlayer(Z),1);    
end
toc











% %%  Trouble-shooting: Draw lines between connected dots
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
%    for z=1:Zs
%             text(DotData2{z}(:,1),DotData2{z}(:,2),[num2str(z)],'color','w','FontSize',8);
%            % text(DotData1{z}(:,1),DotData1{z}(:,2),[num2str(z)],'color','m','FontSize',8);
%    end
%         
% 
% % % Flag specific points of interest.  
% % p1 = (1089+1)/2; p2 = (525);
% % figure(2); hold on; plot(dotC(p1,1),dotC(p1,2),'mo','MarkerSize',15); 
% % plot(dotC(p2,1),dotC(p2,2),'m*','MarkerSize',15);
%    
   
   
   
   
   
   
   
   
   
   
   
   
   

%%

%clear DotMasks DotData R1 Rz Rs LoZ Inds Iw;

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
Cent = cell(1,Nsects); 

for k=1:Nsects
% For each dot in the system, there is a row i, which contains all the dots
% above and below it.  We want to concentrate only on the dots in z~=i
% which are connected to the dot as z=i.  For this we use the LayerJoin
% object.  
     mask = DotConn(1+(k-1)*stp:min(k*stp,2*NDots),:)>0; % where dots exist 
     ConnInt_T = immultiply(ConnInt(1+(k-1)*stp:min(k*stp,2*NDots),:),eval([intype,'(mask)'])); % report intensities where dot connections are found 
     
     % Keep only main diagnol:
     MD = imadd(eval([intype,'(LayerJoin(1+(k-1)*stp:min(k*stp,2*NDots),:))']),ConnInt_T); 
     % figure(4); clf; imagesc(MD);
     MD = MD>0;   
     MD = bwareaopen(MD,20); % Remove all dots in z not connected to the dot at z=i  
     % figure(4); clf; imagesc(MD);
     ConnInt_T = immultiply(ConnInt_T,eval([intype,'(MD)']));
   
    mask = ConnInt_T>0;

%     % Watershed to split dots
%     % W = ConnInt_T.*uint16(mask.*longDots); figure(4); clf; imagesc(W);
%     W = immultiply(ConnInt_T,eval([intype,'(mask)']) ); 
%     % figure(4); clf; imagesc(W);
%     W = watershed(max(W(:)) - W);   % This is the Memory Kill Step; 
%      % figure(3); clf; imagesc(W); colormap lines;
%     mask(W==0) = 0; 
%     % figure(4); clf; imagesc(mask);

    if getpreciseZ == 1
        labeled = bwlabel(mask);
        R = regionprops(labeled,ConnInt_T,'WeightedCentroid');
        tcent =  reshape([R.WeightedCentroid],2,length(R))';

      %  cent = [cent; tcent(:,1),(k-1)*stp + tcent(:,2)];
        Cent{k} = [tcent(:,1),(k-1)*stp + tcent(:,2)];
        
        if plotdata == 1;
            figure(4); clf; 
            colordef black; set(gcf,'color','k'); 
            imagesc( ConnInt_T ); colormap hot; shading flat;  colorbar;
            ylabel('mRNA index'); xlabel('z-depth'); 
            hold on; plot(tcent(:,1),tcent(:,2),'co'); 
            title('Cross-Section of all dots'); 
        end
    end

    % bwareaopen is very expensive for big images
       mask = bwareaopen(mask,consec_layers); % figure(3); clf; imagesc(mask);
       masked_inds(1+(k-1)*stp:min(k*stp,2*NDots),:) =  mask.*DotConn(1+(k-1)*stp:min(k*stp,2*NDots),:)  ;
%      
       
       if comp_onVoff == 1 % For comparison of Dot intensities
        masked_ints(1+(k-1)*stp:min(k*stp,2*NDots),:) =  mask.*single(ConnInt(1+(k-1)*stp:min(k*stp,2*NDots),:) );
       end
end

if plotdata == 1
    cent = cell2mat(Cent); 
end
toc


  %% Troubleshooting
% tic
% 
% figure(2); clf; 
% imagesc(Imax_dots(1:80,1:80,:));   hold on;
% plot(dotC(:,1),dotC(:,2),'c+');
% for n=1:2:2*NDots  
%      text(dotC((n+1)/2,1)+1,dotC((n+1)/2,2),['   ', num2str( (n+1)/2 )],'color','c','FontSize',8);
%     linx = dotC(nonzeros(masked_inds(n,:)),1);
%     liny = dotC(nonzeros(masked_inds(n,:)),2) ;  
%     plot(linx,liny,'w');
% end
% 
% hold on;
%    for z=1:Zs
%             text(DotData2{z}(:,1),DotData2{z}(:,2),[num2str(z)],'color','w','FontSize',8);
%            % text(DotData1{z}(:,1),DotData1{z}(:,2),[num2str(z)],'color','m','FontSize',8);
%    end
%         
% 
% 
% toc









%%
%clear  ConnInt_T DotConn LayerJoin  mask  ConnInt

%%

tic
disp('Building spheres from cross-section disks...');

% figure(3); clf; imagesc(masked_inds); colormap hot;

% testL = NaN*zeros(1,2*NDots); 
remove_dot = zeros(NDots,1); 
stacked_dots =0;
% loop through all dots

for i = 1:2:2*NDots-1 % i = 605 i = 5401; i=5693 i = 5547  i = 6549 
    j = find(masked_inds(i,:));
    counted = masked_inds(i,j(2:end));   
    if isempty(j) == 0
        stacked_dots = max(j)-min(j) > length(j)-1;
    
        if stacked_dots == 0  && getpreciseZ == 1
             ii = find(cent(:,2)==i,1);
             dotC((i+1)/2,3) = cent(ii(1),1);
        end
    else
        remove_dot((i+1)/2) = 1;
    end
    if stacked_dots == 1% if stacked dots split up.  
        brk_pts =[0, find(diff(j)>1),length(j),length(j)]; % breakpoints in stack 
        possibles = masked_inds(i,j); % all possible multicounted indices 
        ii = find(possibles == (i+1)/2,1) ; % find this breakpoint    
        % only need this if low intensity points have been removed
        if isempty(ii); [jnk, ii] = min( ((i+1)/2 - possibles).^2 );  end
       % find nearest breakpoint without going over
          kk = (ii-brk_pts); kk(kk<0) = 100; [jnk,bi] = min(kk);    
          counted = possibles(brk_pts(bi)+2:brk_pts(bi+1)); 
          if  getpreciseZ == 1
              ii = find(cent(:,2)==i);
              dotC((i+1)/2,3) = cent(ii( min(bi,length(ii)) ),1);
          end
    stacked_dots =0;      
    end
     remove_dot(counted) = 1;   
   % testL(i) = length(remove_dot);
end
toc
sum(remove_dot);

% clear DotData DotMasks 
% clear ConnInt ConnInt_T DotConn  LayerJoin R1 Rz Rs  
% clear W Nuclabeld  Cell_bnd dL bw1 
% clear dotC New_dotC test_dotC  Cents
% clear Iin_z LoZ MD conn_map mask masked_inds 

% figure(1); clf; plot(testL,'w.');
% sum(stacked_dots)
          % NB can't sum stacked dots, all stack dots are also multiply
          % counted.  i.e. the first time a doublet is enountered we say
          % Tstacked = Tstacked + 1, and then we enounter that the other of
          % the pair and again say Tstacked = Tstacked + 1;  

 if plotdata ==1
    test_dotC = dotC; test_dotC(logical(remove_dot),3) = 0;
    figure(5); clf; imagesc(masked_inds); colormap hot;
    hold on; plot(test_dotC(:,3),linspace(1,2*NDots-1,NDots),'co');
end
          
          
          
%New_dotC = dotC(~remove_dot,:);
dotC = dotC(~remove_dot,:);
N_dots = length(dotC); 
%N_dots = NDots - sum(remove_dot) % sum(stacked_dots)
%N_dots = length(New_dotC); 
disp(['Counted ',num2str(N_dots),' spheres']); 



if comp_onVoff == 1
    figure(7); clf;  subplot(2,1,1);
    hist(log(nonzeros(masked_ints(~remove_dot,:)))); title('intensities of kept dots');
    subplot(2,1,2);
    hist(log(nonzeros(masked_ints(logical(remove_dot),:)))); title('intensities of removed dots');
end



% figure(2); clf; 
% imagesc(Imax_dots);   hold on;
% plot(dotC(:,1),dotC(:,2),'c+');



