%%  fxnproject.m
% 
%  Alistair Boettiger                                  Date Begun: 11/25/08
% Levine Lab                                        Last Modified: 03/04/10
%
% function script. 
%
%%  Descrition 
% create z-projections from Leica generated tif stacks
% 
%
%% Input variables
% fin       folder containing images
% name      file root
% chns      send layer to chosen color, [R G B].
% subset    project only subseet of the stack, [start_slice, end_slice]
% imnum     image number
% nuc_dat   use different subset for nuclei [nuc_start, nuc_end, nuc_chn] 
%
%%  Output variables
%  Io       output image, 3 layer  
%
%% Called by programs:
% imviewer.m
%
%% Updates
% 
% Updated 03/04/10 to acount for arbitrary 
% Updated on 07/16/10 to change for new Leica format uses z00 not z000. 
% Updated on 08/23/10 to revert to z000



%% Active script

% switch to vargin and use nars 

function Io = fxnproject(fin,name,chns,subset,imnum,nuc_dat)


%% 
% fin = '/Users/alistair/Lab Data/Raw_Data/07.15.09/sogPxdl_25C';
% name = 'sogPxdl_25C_LacZ_sog';
% chns = ['0','1','2'];
% subset = [1,10];
% imnum = '01';
% nuc_dat = [1,10,3];
%%
nuc_start = nuc_dat(1); nuc_end = nuc_dat(2); nuc_chn = nuc_dat(3); 

disp(['fxnproject: ' imnum]); % save test2;

watcher = 0; clear Im;
if subset(2) == 0
    frames = 99;
else
    frames = subset(2) - subset(1)+1; 
end

Im = cell(length(chns),frames);
    
    c = 0; % counter for actual number of frames
    
    for f = subset(1):frames     
        j= subset(1)-1+f;   % 100 is max number of slices used in a projection 

         file = [fin,'/',name,'_',imnum,'_z',num2str(j),'.tif'];
         % save test2;
           try  % try in case the file does not exist
              Imtemp = imreadfast(file); %  imread
            catch er   % stop iterating if exhausted all files.
                disp(er.message);
                watcher = 1;
            end
         
        for i=1:length(chns)
            Im{i}{f} = Imtemp(:,:,i);            
        end
          if watcher == 1; break; end;           
           c=c+1;
    end
    % save test2;
    
   % imsize1 = 1024; imsize2 = 1024;
    [imsize1, imsize2] = size(Im{1}{1}); % get image size
    inttype = class(Im{1}{1});
    Io = zeros(imsize1,imsize2,3,inttype);  % blank to store data
    % convert into stack and max project stack
 
    
    if nuc_end == 0
        nuc_end = c;
    end
      Nchns = length(chns);
       
    %   save test2; 
     
   for i=1:Nchns
       if i ~= nuc_chn
    Itemp = reshape(cell2mat(Im{i}),imsize1,imsize2,length(Im{i}));
     Io(:,:,i) = max(Itemp,[],3);    % mean(Itemp,3); 
       end
   end
 
  nucset = nuc_start:nuc_end; 
  i=nuc_chn;
 % save test2; 
    
    Io(:,:,i) =max(reshape( cell2mat(Im{i}(nucset)),...
        imsize1,imsize2,length(Im{i}(nucset)) ),[],3);
    
    
   figure(21); clf; imshow(Io);   
