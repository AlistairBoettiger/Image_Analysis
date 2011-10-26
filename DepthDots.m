%%                           DepthDots
% Alistair Boettiger                                   Date Begun: 02/07/11        
% Levine  Lab                                       Last Modified: 02/07/11
% 
%% Description
% Written for im_singlemolecule.fig
%
%% Inputs
% In        Image of nuclei (downscaled for speed)
% Cell_bnd  bw image of cell boundaries
% inds_Z    indices of pixels at centroids of localized mRNA transcripts
% h,w       height and width of original image

function DepthDots(In,Cell_bnd,inds_Z,h,w)
                   
% depth color coding of mRNA transcripts
    hn = size(In,1); 
    Zs = length(inds_Z); 

    In = imresize(In,h/hn,'nearest'); % resize nuclei up to image size.
    In = uint8(double(In)/2^16*255); % convert to uint8
    Cell_bnd = uint8(255*imresize(Cell_bnd,h/hn,'nearest'));% resize boundary labels                

    % Add nuclei and cell boundaries to the image
    Idot = cell(1,Zs); 
    Ib = uint8(zeros(h,w,3));
    Ib(:,:,1) = Cell_bnd; 
    Ib(:,:,2) = Cell_bnd; 
    Ib(:,:,3) = In + Cell_bnd;
    Iv = Ib;
    col = spring(Zs); % colormap for depth coding
   for z=1:Zs;
       I1 = false(h,w);
       I1(inds_Z{z}) = 1; % place all dots on array
       % Paint different color for dots of each z-plane.  
       Iv(:,:,1) = Iv(:,:,1) + uint8(col(z,1)*I1*255);  
       Iv(:,:,2) =Iv(:,:,2) +  uint8(col(z,2)*I1*255);
       Iv(:,:,3) = Iv(:,:,3) + uint8(col(z,3)*I1*255);
       Idot{z} = I1; 
   end
     figure(7); clf; colormap(col);
     colordef black; set(gcf,'color','k'); 
     imshow(Iv); colorbar; caxis([1,Zs]);
     
     