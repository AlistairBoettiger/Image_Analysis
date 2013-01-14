% cyto2nucmap
% use cytoplasmically distributed mRNA to extract nuclei map.  


folder = 'D:\Data\2013-01-10_no-dist\';
fout = folder;
fname = 'MP08Hz_snaD_2';

nuc_chn = 3; 


    %%
    emb = 6;

C = zeros(2048,2048,'uint16');

zs = 30; % starting frame for z-projection of nuclei channel
zf = 35;  % ending frame fo z-projection of nuclei channel

for z =   zs:zf;
    I = imreadfast([folder,fname, '_0', num2str(emb),'_z',num2str(z),'.tif']); 
   % figure(1); clf; imshow(I);
    C = (C+I(:,:,nuc_chn))./uint16(2);
end

F = fspecial('gaussian',50,5); 

C2 = imfilter(C,F,'replicate');
C2 = imresize(C2,512/2048);
figure(1); clf; imagesc(C2); colormap gray; colorbar;
%%

N = imadjust(C2,[0,1],[0,1],.45);
 N = max(N(:)) - N;

% N = .1*max(C2(:)) - C2;
figure(2); clf; imagesc(N); colormap gray; colorbar; caxis([0,1*max(N(:))])

%N = C2;
%
% minN = 100; imblur = 2; sigmaE = 26; sigmaI = 30; FiltSize = 45;
% Mthink = 45; Mthin = 3; Imnn = 2;
minN = 40; imblur = 2; sigmaE = 15; sigmaI = 17; FiltSize = 45;
Mthink = 45; Mthin = 3; Imnn = 2;
 [handles.bw,handles.cent] = fxn_nuc_seg(N,minN,sigmaE,sigmaI,FiltSize,imblur);
%%
 [NucLabeled,Nuc_overlay,conn_map,cell_bords] = fxn_nuc_reg(N,handles.bw,Mthink,Mthin,Imnn);  
  
  [h,w] = size(NucLabeled);
    
    Cell_bnd = false(h,w);
    Cell_bnd(cell_bords) = 1;
    Cell_bnd = bwareaopen( Cell_bnd,100);
    
     handles.Nucs = imresize(I(:,:,1),512/2048);
    
    int = class(handles.Nucs);
    if strcmp(int,'uint16')
        Nucs = handles.Nucs + uint16(2^16*Cell_bnd);
    elseif strcmp(int,'uint8')
         Nucs = handles.Nucs + uint8(2^8*Cell_bnd);
    end
        
    
    figure(21); clf; imagesc(Nucs); colormap('hot');   
   %% 
    nuc_cents = handles.cent; 
    Nucs = N; 
    AgeClass = input('embryo age class =','s');
    
    save([fout,'/',fname,'_0',num2str(emb),'_nucdata.mat'],...
        'NucLabeled','nuc_cents','conn_map','Cell_bnd','Nucs','AgeClass'); 
    %%
