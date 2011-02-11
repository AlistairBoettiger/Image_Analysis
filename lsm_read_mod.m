
% Alistair Boettiger
% Levine lab
% UC Berkeley 


function stack = lsm_read_mod(name,N)

%% 
%  clear all;
%  name = '/Volumes/Data/Lab Data/Raw_Data/02-08-11/MP10_22C_sna_y.mat';
%  N =9; 

global TIF;
load(name)
filetemp= fopen(Datas.filename,'r','l');
Zs = Datas.LSM_info.DimensionZ;
filetemp= fopen(Datas.filename,'r','l');

stack= cell(1,Zs);

for i=1:Zs  % i = 21

    TIF = Datas.([ 'Stack' num2str(N)]).(['Image' num2str(i)]).TIF;
    IMG = Datas.([ 'Stack' num2str(N)]).(['Image' num2str(i)]).IMG;
    TIF.file=filetemp;


    
           offset = 0; 
            %read the image channels
            for c = 1:TIF.SamplesPerPixel % c=3
                TIF.StripCnt = c;
                IMG.data{c} = read_planeT(offset, IMG.width, IMG.height, c,TIF); 
          
                % check for screwed up offset
                [h,w] = size(IMG.data{c});
                % look in middle of data set and see if it's noise or
                % signal
                sdata =  IMG.data{c}( floor(h/2*.9):floor(h/2*1.1), floor(w/2*.9):floor(w/2*1.1)  ); 
                isnoise = std(double(sdata(:)));
               if isnoise> 1.4E4 
                   offset = 1;
                   IMG.data{c} = read_planeT(offset, IMG.width, IMG.height,c,TIF); 
                   sdata =  IMG.data{c}( floor(h/2*.9):floor(h/2*1.1), floor(w/2*.9):floor(w/2*1.1)  ); 
                   fixed = std(double(sdata(:)));
                   disp(['offset error found in chn ', num2str(c), ' layer ',num2str(i),...
                       '  std=',num2str(isnoise,5), '  now=',num2str(fixed,5)] );
                   
                   if fixed > 1.4E4
                       if TIF.BitsPerSample(1) == 16
                          IMG.data{c} = uint16(zeros(h,w));
                       else
                          IMG.data{c} = uint8(zeros(h,w));
                       end
                   end
                  disp(['Fix failed for in chn ', num2str(c), ' layer ',num2str(i),...
                       '  skipping this image...'] );
                   
              end
           end           
            stack{1,i} = IMG.data; 
       
end
 

% Just for troubleshooting
% 
% for z = 1:Zs
%     figure(1); clf;  imshow(stack{1,21}{3});
%     pause(.01);
% end

% clear all;

