% cyto2nucmap
% use cytoplasmically distributed mRNA to extract nuclei map.  


folder = 'C:\Users\Alistair\Data\2011-12\';
fout = folder;
fname = 's05_MP08_b';

for emb = 5;
    %%
    emb = 2;

C = zeros(2048,2048,'uint16');

for z =   22:30;
    I = imread([folder,fname, '_0', num2str(emb),'_z',num2str(z),'.tif']); 
    C = (C+I(:,:,2))./uint16(2);
end

F = fspecial('gaussian',50,5); 

C2 = imfilter(C,F,'replicate');
C2 = imresize(C2,512/2048);
figure(1); clf; imagesc(C2); colormap gray; colorbar;
%%

N = imadjust(C2,[0,.061],[0,1],.51); N = max(N(:)) - N;

% N = .1*max(C2(:)) - C2;
figure(2); clf; imagesc(N); colormap gray; colorbar; caxis([0,1*max(N(:))])

%N = C2;

minN = 100; imblur = 2; sigmaE = 16; sigmaI = 18; FiltSize = 30;
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
    
    save([fout,'/',fname,'_0',num2str(emb),'_nucdata.mat'],...
        'NucLabeled','nuc_cents','conn_map','Cell_bnd','Nucs'); 
    %%
end