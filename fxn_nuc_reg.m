%%                                  fxn_nuc_reg.m                        %%
% Alistair Boettiger & Jacques Bothma           Version Completed: 01/28/10          
% Levine Lab                                        Last Modified: 01/28/10

% This picks up from fxn_nuc_seg, which segements the nuclei

function  [Nuc_map,Nuc_cents,conn_map,A] = fxn_nuc_reg(I,bw,Mthink,Mthin,Imnn)

% % Trouble shooting defaults.  
%  clear all; close all;
%  
%  % % % save nuc_seg_testdata I bw;
%  load nuc_seg_testdata; 
% Mthink = 15;
% Mthin = 3;
% Imnn = 2; 
%      

%% Fill Gaps
% this section uses Mthink, Mthin, Imnn

          % for save data; 
fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';
          %%%%%%%% code that splits the nuclei of the embryo into seperate segments
          
L = bwlabeln(bw,8);

EA = imfill(imopen(imdilate(I>0,strel('disk',2)),strel('disk',5)),'holes'); %EA is a logical array that is 1 in the rough region that corresponds to the embryo

LL = bwmorph(L,'thicken',Mthink); % This function expands the individual nuclei by adding pixels to the outer edge untill the nuclei come within 1 pixel
                                  % of touching or M think pixels are
                                  % added.
% save([fdata,'/','test2']);
LL=LL.*EA;                        % Line that trims the extra expansion beyond embryo boundary

EAN=imclose(LL,strel('disk',Imnn)); 
EAN=imfill(EAN,'holes');             %Better definition of embryo region.


K = bwmorph(~LL.*EAN,'thin',Mthin);   % Code that thins the lines sepearting individual nuclear segments to a single pixel.

[La,NUM] = bwlabeln(EAN-K,4);          %Code that labels unique nuclear segments 


% figure, imshow(label2rgb(La, 'jet', [1,1,1],'shuffle'),'Border','tight','InitialMagnification',100);% maxwindow

LaB=La>0;

[h,w] = size(K);  

A=find(K==1);  % figure(3); clf; imshow(K); 


%%%%%%Code that assigns the edge pixels seperating nuclear segments to distinct nuclear
%%%%%%regions based on neighbouring pixel intesities. Note the use of linear indexing which speeds things up. 


AA=size(I,1);

A = setdiff(A,1:AA); 
A = setdiff(A,h*w-AA:h*w); 
N = length(A); 


MV=[AA,AA+1,1,-AA+1,-AA,-AA-1,-1,AA-1]; %8 closest neighbours
%MV2=[-2*AA-2,-AA-2,-2,AA-2,2*AA-2,-2*AA+2,-AA+2,+2,AA+2,2*AA+2,2*AA-1,2*AA,2*AA-1]; %16 next nearest neighbours
%MV=[MV,MV2];
La2=La;

for i=1:N
    sindex=A(i);   
    idx= MV+sindex;
    
    PV1=double(I(idx)');  % intensities of neighbor pixels
    PV2=nonzeros(La(idx)');  % labels of neighbor pixels 
    try
        La2(sindex)=PV2(end);
    catch
        'empty';  
    end
end


% In order to figure out which nuclear domains are next to one another you
% need to determine which of the elements in your label array is within one
% pixel of another element that is different to it. Instead of looping
% through all elements in the label matrix I circularly shift the whole
% label matrix by one (again you need to have a buffer around the edge not
% to run into strife) pixel in each direction (That is what S is for) and
% then take the difference between the original label matrix and the
% shifted one.  The only elements that are non-zero in this matrix are the
% ones where the label of the element is different to the pixel next to it
% in the direction in which the label matrix had been shifted. These two
% elements then constitute a place where two different nuclear domains are
% next to one another and so by simply finding the labels of these two 
% elements (linear index of the shifted on stored in MV) you can create a
% two element row vector that corresponds to nuclei that are touching.
% Putting all of these into a big column vector gives all the
% connectivities through boundaries and you convert that into a
% conventional connectivity matrix using the accumarray command. It not
% only keeps track of which regions are connected but also through how many
% elements. Will specify exactly what's going on line by line in red. Let
% me know if you want further clarification.





    %%%%%%%%% Defining connectivity of nuclear domains %%%%%%%%

    S=[ 0 1; 1 0; -1 0; 0 -1; 1 1; -1 1; 1 -1; -1 -1]; %Array that specifies the 8 possible directions to move the array in to cover all nearest neighbours.
    AA=size(La2,1);  %Specifies the number of rows of the label matric that will be used in the linear indexing

    MV=[AA,1,-1,-AA,AA+1,+AA-1,-AA+1,-AA-1];  %Specifies the numbers that need to be added to the linear index of an element to determine the lionear index of the neighbour that is reached when the matrix is

                                                                                     %circularly shifted by the correponding element in S
    NA=[];  NB=[];  
    La2M=La2;    La2M(La2M==0)=NaN; % putting nans in the label matrix where there is zeros to make sure nuclei next to the background isn't linked to it.
[h,w] = size(La2M);
 La2M(1:end,1) = NaN; La2M(1:end,end) = NaN; 
 La2M(1,1:end) = NaN; La2M(end,1:end) = NaN; 

    for i=1:8  % Loop that goes through all nearest neighbour pairs
       NAT=find(abs((La2M-circshift(La2M,S(i,:))))>0); %Taking difference of label matrix and it's shifted version and returning linear indices of non zero elements
       NA=[NA;NAT]; %Store these indices
       NB=[NB;NAT-MV(i)]; %Detemine and store the indices of the neighbouring element
    end

% save([fdata,'/','test2']);
    NAV=La2(NA);  %Label of one member of nuclear pair
    NBV=La2(NB); %Label of other member of nuclear pair
    NN=[NAV,NBV];   %Column vector with all pairs of nuclear linkages
    ACON = accumarray(NN,1); %command that turns these into a connectivity matrix



%%%%%% Display  %%%%%%%%%
                CM=label2rgb(La2, 'jet', [1,1,1],'shuffle');    
               % figure, imshow(CM);

         DI = uint8(bsxfun(@times,double(CM)/255,double(I)));
        figure(21); imshow(DI,'Border','tight','InitialMagnification',100);% maxwindow

 % compute centroids       
        S = regionprops(La2,'Centroid'); % measure areas
cent = reshape([S.Centroid],2,length(S));

% export
Nuc_map = La2; 
Nuc_cents = cent; 
conn_map = ACON;



