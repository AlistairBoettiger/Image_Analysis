
%                       connectivity.m
%
%
% Jacques Bothma and Alistair Boettiger
%
% Description:
% Takes a space filled labeled image and returns the number of connected
% pixels between each object and all other objects.  
%
% Inputs
% L is a labeled matrix
%
% Outputs 
% conn_map   is an N x N matrix (where N is the number of objects in L)
%               where the number of connections between object i and object
%               j is in position i,j.  

function conn_map  = connectivity(L)


%%%%%%%%% Defining connectivity of nuclear domains %%%%%%%%

    S=[ 0 1; 1 0; -1 0; 0 -1; 1 1; -1 1; 1 -1; -1 -1]; %Array that specifies the 8 possible directions to move the array in to cover all nearest neighbours.
    AA=size(L,1);  %Specifies the number of rows of the label matric that will be used in the linear indexing

    MV=[AA,1,-1,-AA,AA+1,+AA-1,-AA+1,-AA-1];  %Specifies the numbers that need to be added to the linear index of an element to determine the lionear index of the neighbour that is reached when the matrix is

                                                                                     %circularly shifted by the correponding element in S
    NA=[];  NB=[];  
    La2M=L;   
    
    La2M(La2M==0)=NaN; % putting nans in the label matrix where there is zeros to make sure nuclei next to the background isn't linked to it.

 La2M(1:end,1) = NaN; La2M(1:end,end) = NaN; 
 La2M(1,1:end) = NaN; La2M(end,1:end) = NaN; 

    for i=1:8  % Loop that goes through all nearest neighbour pairs
       NAT=find(abs((La2M-circshift(La2M,S(i,:))))>0); %Taking difference of label matrix and it's shifted version and returning linear indices of non zero elements
       NA=[NA;NAT]; %Store these indices
       NB=[NB;NAT-MV(i)]; %Detemine and store the indices of the neighbouring element
    end

% save([fdata,'/','test2']);
    NAV=L(NA);  %Label of one member of nuclear pair
    NBV=L(NB); %Label of other member of nuclear pair
    NN=[NAV,NBV];   %Column vector with all pairs of nuclear linkages
    conn_map = accumarray(NN,1); %command that turns these into a connectivity matrix
