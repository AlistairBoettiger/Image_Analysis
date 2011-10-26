%%                          checkZsep
% Alistair Boettiger                                   Date Begun: 02/07/11        
% Levine  Lab                                       Last Modified: 02/07/11
% 
% Written for im_singlemolecule.fig
% Optional plotting command: shows adjacent z-section and flags dots ID'd as
% unique to that section. 

%% Inputs:
%  z          layer to look at
%  Cell_bnd   bw image of determined boundaries of individual cells
%  Im         raw data in mRNA channel. Each cell contains a different z plane  
%  DotData    x,y centroid coordinates of dots localized in plane z
%  D2u_Z      x,y centroids of *unique* dots in layer z.


function checkZsep(z,Cell_bnd,Im,mRNAchn,DotData,D2u_Z)

    [h,w] = size(Im{1,1}{1}); 
    hn = size(Cell_bnd,1);
 
    Cell_bnd = uint16(2^16*imresize(Cell_bnd,h/hn,'nearest'));
      Iz = uint16(zeros(h,w,3));
      Iz(:,:,1) = 2*Im{1,z-1}{mRNAchn} + Cell_bnd;
      Iz(:,:,2) = 2*Im{1,z-2}{mRNAchn} + Cell_bnd;
      Iz(:,:,3) = 2*Im{1,z}{mRNAchn} + Cell_bnd;

       figure(5); clf;  
       imshow(Iz);    hold on;   
       plot(DotData{z-2}(:,1),DotData{z-2}(:,2),'y+');
       plot(DotData{z-1}(:,1),DotData{z-1}(:,2),'ro');
       plot(DotData{z}(:,1),DotData{z}(:,2),'co'); 
       plot(D2u_Z{z-1}(:,1),D2u_Z{z-1}(:,2),'r.','MarkerSize',10);
       plot(D2u_Z{z}(:,1),D2u_Z{z}(:,2),'c.','MarkerSize',10);
