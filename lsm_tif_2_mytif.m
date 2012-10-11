
% convert Zeiss Tifs labels to my Tifs

clear all;

folder = 'D:\Data\2012-08-20_ML124\' ; % 'D:\Data\2012-08-15_MLslides\';%   %  'D:\Data\2012-08-26_mp07\'; % 
subfolder =  'm124\'; % 'm106a\';  % 'mp07b\'; %   
froot ='ml124B'; % 'm106a'; %  'mp07b'; % 

% dir(folder); 

Zs =65;
Ims = 10;
r = 2; % channel to make red
g = 3; % channel to make green; 
b = 1; % channel to make blue



for i=10:Ims
    Imax = zeros(2048,2048,3,'uint16'); % so much for flexible image size
    for z=1:Zs
        fname = [froot,'_z',sprintf('%02d',z-1),'_p',sprintf('%02d',i-1),'.tif'];
        try
            Iin = imreadfast([folder,subfolder,fname]);
        catch er
            disp(er.message)
            continue
        end
        [h,w,chns] = size(Iin);
        intype = class(Iin); 
        Iout = zeros(h,w,3,intype);
        if r<= chns
            Iout(:,:,1) = Iin(:,:,r); 
            Imax(:,:,1) = max( cat(3,Imax(:,:,1),Iin(:,:,r)),[],3);   
        end
        if g<= chns
            Iout(:,:,2) = Iin(:,:,g); 
            Imax(:,:,2) = max( cat(3,Imax(:,:,2),Iin(:,:,g)),[],3);   
        end
        if b <= chns
            Iout(:,:,3) = Iin(:,:,b); 
            Imax(:,:,3) = max( cat(3,Imax(:,:,3),Iin(:,:,b)),[],3);   
        end
        figure(1); clf; 
        subplot(1,2,1); imagesc(Iout);
        subplot(1,2,2); imagesc(Imax); 
        pause(.1);
        
        fout = [froot,'_',sprintf('%02d',i),'_z',num2str(z),'.tif'];
        imwrite(Iout,[folder,fout]); 
        disp([num2str(i),': ',num2str(z/Zs*100,2),'% finished']);
    end
        imax_name = ['max_',froot,'_',sprintf('%02d',i),'.tif'];
    imwrite(Imax,[folder,imax_name]); 
end
