

output_folder = 'C:\Users\Alistair\Documents\Projects\Chromatin\Results';
input_folder = 'D:\Data\2012-11-14_MTDs';
% fname = 'max_SCM_late_CNS_01';
% fname = 'max_SCM_late_CNS_ectoderm_01';
fname = 'max_SCM_onsetGBE2_abdB_iab8_K27_01';

I = imread([input_folder,filesep,fname,'.tif']);

figure(1); clf; Ncolor(I);


cmap = [.25,.25,.25;
        1, 0, 0;  % iab8
        0,.7, .7; % abdA
        .7,0,.7];
    
    figure(1); clf; Ncolor(I,cmap);
    
    I(:,:,1) = mycontrast(I(:,:,1),.001,.001);
    I(:,:,2) = mycontrast(I(:,:,2),.00001,.95);
    I(:,:,3) = mycontrast(I(:,:,3),.0001,.001);
  %  I(:,:,3) = mycontrast(I(:,:,3),.00001,.25);
 % I(:,:,4) = mycontrast(I(:,:,4),.0001,.0001);
     Iout = Ncolor(I,cmap);
    figure(1); clf; imagesc(Iout);
    
    % save smaller tifs of image 
    imwrite(imresize(Iout,.2),[output_folder,filesep,fname,'.png']);
 imwrite(imresize(I(:,:,1),.2),[output_folder,filesep,fname,'_c1.png']);
 imwrite(imresize(I(:,:,2),.2),[output_folder,filesep,fname,'_c2.png']);
 imwrite(imresize(I(:,:,3),.2),[output_folder,filesep,fname,'_c3.png']);
 imwrite(imresize(I(:,:,4),.2),[output_folder,filesep,fname,'_c4.png']);
