
% Alistair Boettiger
% Levine lab
% UC Berkeley 


function stack = lsm_read_mod(name,N)

% % 
%  clear all;
%  name = '/Volumes/Data/Lab Data/Raw_Data/02-08-11/MP10_22C_sna_y.mat';
%  N =5; i = 11; 

global TIF;
load(name)
filetemp= fopen(Datas.filename,'r','l');
Zs = Datas.LSM_info.DimensionZ;
filetemp= fopen(Datas.filename,'r','l');

stack= cell(1,Zs);

for i=1:Zs

    TIF = Datas.([ 'Stack' num2str(N)]).(['Image' num2str(i)]).TIF;
    IMG = Datas.([ 'Stack' num2str(N)]).(['Image' num2str(i)]).IMG;
    TIF.file=filetemp;

           offset = 0; 
            %read the image channels
            for c = 1:TIF.SamplesPerPixel
                TIF.StripCnt = c;
                IMG.data{c} = read_planeT(offset, IMG.width, IMG.height, c,TIF); 
          
                % check for screwed up offset
                sdata =  IMG.data{c}(1:100,1:100); 
                isnoise = std(double(sdata(:)));
               if isnoise> 5E3 
                   offset = 1;
                   IMG.data{c} = read_planeT(offset, IMG.width, IMG.height,c,TIF); 
               end
           end           
            stack{1,i} = IMG.data; 
       
end
 

% % Just for troubleshooting
%   figure(1); clf; 
% for z = 1:Zs
%   imshow(stack{1,z}{1});
%     pause(.01);
% end
% 
