%% mRNA dynamics, show different nuclei patterns. 
% Alistair Boettiger                                Date Begun: 12/20/11
% Zhuang Lab / Levine Lab                           Last Modified: 07/31/12


%% Publication figure as of 09/26/12

%%
% load all data sorted by age
clear all;

   folder = 'C:\Users/Alistair/My Documents/Projects/mRNA_counting/Data/';
   rawfolder = 'D:/Data/';
   rawfolder1 = 'G:/Raw_Data/'; 
   
   savedata = 0; % don't write data to disk (will auto-overwrite existing)
   loadnucdata = 0;
%folder = '/home/alistair/Documents/Research/Projects/mRNA_counting/Data/';
% rawfolder = '/media/ALISTAIR2/Data/';
 %% Load Nuc data  
 for n=1:1
     Nuc = cell(57,1); 
 if loadnucdata == 1
 % ----------------------------% exiting cc12 ---------------------------%
% 'wt_sna'  emb = 03, 04, 10
% 'MP06_cflip' emb = 08
% 'MP06Hz_b' emb = 02

Nuc{5} = imread([rawfolder,slidedate,'max_',fname,'_03.tif']);
Nuc{3} = imread([rawfolder,'2011-12/','max_wt_sna_04.tif']);
Nuc{1} = imread([rawfolder,'2011-12/','max_wt_sna_10.tif']);
Nuc{4} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_08.tif']);
Nuc{2}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_02.tif']);
%-----------------------------------------------------------------------%

%========================== during cc13====================================
% 'wt_sna'  emb = 09
% 'MP06_cflip' emb = 01, 02, 09, 10 
% 'MP06Hz_b' emb = 04, 05

Nuc{6} = imread([rawfolder,'2011-12/','max_wt_sna_09.tif']);
Nuc{7} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_01.tif']);
Nuc{8} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_02.tif']);
Nuc{12} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_09.tif']);
Nuc{11} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_10.tif']);
Nuc{9}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_04.tif']);
Nuc{10}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_05.tif']);
%========================================================================%


% --------------------- % exiting cc13 -------------------------------------
% 'wt_sna'  emb = 05, 06 
% 'MP06_cflip' emb = 03, 04, 06
% 'MP06Hz_b' emb = 01, 09, 10
Nuc{15} = imread([rawfolder,'2011-12/','max_wt_sna_05.tif']);  % meta
Nuc{14} = imread([rawfolder,'2011-12/','max_wt_sna_06.tif']); % meta
Nuc{13} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_03.tif']); % meta
Nuc{20} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_04.tif']); % ana
Nuc{17} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_06.tif']); %meta
Nuc{18}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_01.tif']); % late meta
Nuc{19}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_09.tif']); % late meta
Nuc{16}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_10.tif']); % ana
%------------------------------------------------------------------%

% =============== % telophase into cc14 ===============================
% 'MP06_cflip' emb = 05 
% 'MP06Hz_b' emb = 03, 07
% 'MP06Hz' emb = 10
Nuc{21}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_07.tif']);
Nuc{22}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_10.tif']);
Nuc{23} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_05.tif']);
Nuc{24}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_03.tif']);
% =====================================================================



% % ---------------- early cc14 ---------------------------%
% 'wt_sna'  emb = 01, 02,  07, 08,
% 'MP06Hz_b' emb = 06, 08, 11
% 'MP06Hz'; emb = 01, 02, 11, 12
Nuc{25} = imread([rawfolder,'2011-12/','max_wt_sna_01.tif']);
Nuc{26} = imread([rawfolder,'2011-12/','max_wt_sna_02.tif']);
Nuc{27} = imread([rawfolder,'2011-12/','max_wt_sna_07.tif']);
Nuc{32} = imread([rawfolder,'2011-12/','max_wt_sna_08.tif']);
Nuc{29}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_06.tif']);
Nuc{30}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_08.tif']);
Nuc{31}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_11.tif']);
Nuc{28}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_01.tif']);
Nuc{33}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_02.tif']);
Nuc{34}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_11.tif']);
Nuc{35}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_12.tif']);
%---------------------------------------------------------%


% ========================  steady state cc14 ~170 ==============================% 
% 'MP06Hz' emb = 03, 04, 05, 06, 07, 08
Nuc{36}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_07.tif']);
Nuc{37}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_06.tif']);
Nuc{56}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_03.tif']);
Nuc{57}= imread([rawfolder1,'2011-05-22/s05_MP06/','max_s05_MP06Hz_04.tif']);

% =================================================================== % 
 end

 end

   
 %% Load mRNA data  
for n=1:1 
 clear mRNAs mu; 
% ----------------------------% exiting cc12 ---------------------------%
% 'wt_sna'  emb = 03, 04, 10
% 'MP06_cflip' emb = 08
% 'MP06Hz_b' emb = 02
slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{5} = data{3}.Data_sort(:,2);  mu12telo{5} = data{3}.mu;
mRNAs12telo{3} = data{4}.Data_sort(:,2);  mu12telo{3} = data{4}.mu;
x12telo{5} = data{3}.Data_sort(:,1)/1000;
x12telo{3} = data{4}.Data_sort(:,1)/1000;
% mRNAs12telo{1} = data{10}.Data_sort(:,2); mu12telo{1} = data{10}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{4} = data{8}.Data_sort(:,2); mu12telo{4} = data{8}.mu;
x12telo{4} = data{8}.Data_sort(:,1)/1000;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{2} = data{2}.Data_sort(:,2); mu12telo{2} = data{2}.mu;
x12telo{2} = data{2}.Data_sort(:,1)/1000;

slidedate = '2011-12/'; fname = 's142_sna';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{6} = data{1}.Data_sort(:,2); mu12telo{6} = data{1}.mu;
mRNAs12telo{7} = data{2}.Data_sort(:,2); mu12telo{7} = data{2}.mu;
x12telo{6} = data{1}.Data_sort(:,1)/1000;
x12telo{7} = data{2}.Data_sort(:,1)/1000;
%-----------------------------------------------------------------------%

%========================== during cc13====================================
% 'wt_sna'  emb = 09
% 'MP06_cflip' emb = 01, 02, 09, 10 
% 'MP06Hz_b' emb = 04, 05
slidedate = '2011-12/'; fname = 's142_sna';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{5} = data{10}.Data_sort(:,2); mu13{5} = data{10}.mu;
mRNAs13{9} = data{11}.Data_sort(:,2); mu13{9} = data{11}.mu;
mRNAs13{1} = data{8}.Data_sort(:,2); mu13{1} = data{8}.mu;
x13{5} = data{10}.Data_sort(:,1)/1000;
x13{9} = data{11}.Data_sort(:,1)/1000;
x13{1} = data{8}.Data_sort(:,1)/1000;

slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{2} = data{9}.Data_sort(:,2); mu13{2} = data{9}.mu;
x13{2} = data{9}.Data_sort(:,1)/1000;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{3} = data{1}.Data_sort(:,2); mu13{3} = data{1}.mu;
mRNAs13{4} = data{2}.Data_sort(:,2); mu13{4} = data{2}.mu;
mRNAs13{8} = data{10}.Data_sort(:,2); mu13{8} = data{10}.mu;
mRNAs13{10} = data{9}.Data_sort(:,2); mu13{10} = data{9}.mu;
x13{3} = data{1}.Data_sort(:,1)/1000;
x13{4} = data{2}.Data_sort(:,1)/1000;
x13{8} = data{10}.Data_sort(:,1)/1000;
x13{10} = data{9}.Data_sort(:,1)/1000;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{6} = data{4}.Data_sort(:,2); mu13{6} = data{4}.mu;
mRNAs13{7} = data{5}.Data_sort(:,2); mu13{7} = data{5}.mu;
x13{6} = data{4}.Data_sort(:,1)/1000;
x13{7} = data{5}.Data_sort(:,1)/1000;
%========================================================================%



% --------------------- % metaphase cc13 ----------------------------------
% 'wt_sna'  emb = 05, 06 
% 'MP06_cflip' emb = 03, 04, 06
% 'MP06Hz_b' emb = 01, 09, 10

slidedate = '2011-12/'; fname = 's142_sna';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{1} = data{15}.Data_sort(:,2); mu13meta{1} = data{15}.mu;
x13meta{1} = data{15}.Data_sort(:,1)/1000;

slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{4} = data{5}.Data_sort(:,2); mu13meta{4} = data{5}.mu;
mRNAs13meta{3} = data{6}.Data_sort(:,2); mu13meta{3} = data{6}.mu;
x13meta{4} = data{5}.Data_sort(:,1)/1000;
x13meta{3} = data{6}.Data_sort(:,1)/1000;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{2} = data{3}.Data_sort(:,2); mu13meta{2} = data{3}.mu;
%mRNAs13meta{9} = data{4}.Data_sort(:,2); mu13meta{9} = data{4}.mu;
x13meta{2} = data{3}.Data_sort(:,1)/1000;
%x13meta{9} = data{4}.Data_sort(:,1)/1000;

slidedate = '2011-12/';   fname ='s140_sna'; ver = '_v2';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data');
mRNAs13meta{6} = data{7}.Data_sort(:,2); mu13meta{6} = data{7}.mu;
x13meta{6} = data{7}.Data_sort(:,1)/1000;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{7} = data{1}.Data_sort(:,2); mu13meta{7} = data{1}.mu;
mRNAs13meta{8} = data{9}.Data_sort(:,2); mu13meta{8} = data{9}.mu;
mRNAs13meta{5} = data{10}.Data_sort(:,2); mu13meta{5} = data{10}.mu;
x13meta{7} = data{1}.Data_sort(:,1)/1000;
x13meta{8} = data{9}.Data_sort(:,1)/1000;
x13meta{5} = data{10}.Data_sort(:,1)/1000;
%------------------------------------------------------------------%

% =============== % telophase into cc14 ===============================
% 'MP06_cflip' emb = 05 
% 'MP06Hz_b' emb = 03, 07
% 'MP06Hz' emb = 10
slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{3} = data{5}.Data_sort(:,2); mu13telo{3} = data{5}.mu;
x13telo{3} = data{5}.Data_sort(:,1)/1000;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{4} = data{3}.Data_sort(:,2); mu13telo{4} = data{3}.mu;
mRNAs13telo{1} = data{7}.Data_sort(:,2); mu13telo{1} = data{7}.mu;
x13telo{4} = data{3}.Data_sort(:,1)/1000;
x13telo{1} = data{7}.Data_sort(:,1)/1000;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{2} = data{10}.Data_sort(:,2);  mu13telo{2} = data{10}.mu;
x13telo{2} = data{10}.Data_sort(:,1)/1000;
% =====================================================================



% % ---------------- early cc14 ---------------------------%
% 'wt_sna'  emb = 01, 02,  07, 08,
% 'MP06Hz_b' emb = 06, 08, 11
% 'MP06Hz'; emb = 01, 02, 11, 12
slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14e{2} = data{1}.Data_sort(:,2);   mu14e{2} = data{1}.mu;
mRNAs14e{4} = data{2}.Data_sort(:,2);   mu14e{4} = data{2}.mu;
mRNAs14e{10} = data{7}.Data_sort(:,2);   mu14e{10} = data{7}.mu;
mRNAs14e{11} = data{8}.Data_sort(:,2);   mu14e{11} = data{8}.mu;
x14e{2} = data{1}.Data_sort(:,1)/1000;
x14e{4} = data{2}.Data_sort(:,1)/1000;
x14e{10} = data{7}.Data_sort(:,1)/1000;
x14e{11} = data{8}.Data_sort(:,1)/1000;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14e{5} = data{6}.Data_sort(:,2);   mu14e{5} = data{6}.mu;
mRNAs14e{14} = data{8}.Data_sort(:,2);  mu14e{14} = data{8}.mu;
mRNAs14e{9} = data{11}.Data_sort(:,2);  mu14e{9} = data{11}.mu;
x14e{5} = data{6}.Data_sort(:,1)/1000;
x14e{14} = data{8}.Data_sort(:,1)/1000;
x14e{9} = data{11}.Data_sort(:,1)/1000;

slidedate = '2011-12/';   fname ='s140_sna'; ver = '_v2';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data');
mRNAs14e{12} = data{1}.Data_sort(:,2); mu14e{12} = data{1}.mu(:,1);
mRNAs14e{3} = data{3}.Data_sort(:,2); mu14e{3} = data{3}.mu(:,1);
mRNAs14e{15} = data{4}.Data_sort(:,2); mu14e{15} = data{4}.mu(:,1);
mRNAs14e{6} = data{5}.Data_sort(:,2); mu14e{6} = data{5}.mu(:,1);
mRNAs14e{13} = data{8}.Data_sort(:,2); mu14e{13} = data{8}.mu(:,1);
mRNAs14e{16} = data{9}.Data_sort(:,2); mu14e{16} = data{9}.mu(:,1);
x14e{12} = data{1}.Data_sort(:,1)/1000;
x14e{3} = data{3}.Data_sort(:,1)/1000;
x14e{15} = data{4}.Data_sort(:,1)/1000;
x14e{6} = data{5}.Data_sort(:,1)/1000;
x14e{13} = data{8}.Data_sort(:,1)/1000;
x14e{16} = data{9}.Data_sort(:,1)/1000;

slidedate = '2011-12/'; fname = 's142_sna';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14e{19} = data{3}.Data_sort(:,2); mu14e{19} = data{3}.mu;
mRNAs14e{18} = data{4}.Data_sort(:,2); mu14e{18} = data{4}.mu;
mRNAs14e{20} = data{5}.Data_sort(:,2); mu14e{20} = data{5}.mu;
mRNAs14e{8} = data{6}.Data_sort(:,2); mu14e{8} = data{6}.mu;
mRNAs14e{7} = data{7}.Data_sort(:,2); mu14e{7} = data{7}.mu;
mRNAs14e{17} = data{9}.Data_sort(:,2); mu14e{17} = data{9}.mu;
mRNAs14e{1} = data{12}.Data_sort(:,2); mu14e{1} = data{12}.mu;
mRNAs14e{21} = data{13}.Data_sort(:,2); mu14e{21} = data{13}.mu;
x14e{19} = data{3}.Data_sort(:,1)/1000;
x14e{18} = data{4}.Data_sort(:,1)/1000;
x14e{20} = data{5}.Data_sort(:,1)/1000;
x14e{8} = data{6}.Data_sort(:,1)/1000;
x14e{7} = data{7}.Data_sort(:,1)/1000;
x14e{17} = data{9}.Data_sort(:,1)/1000;
x14e{1} = data{12}.Data_sort(:,1)/1000;
x14e{21} = data{13}.Data_sort(:,1)/1000;
%---------------------------------------------------------%

% ========================  steady state cc14 ~170 ==============================% 
% 'MP06Hz' emb = 03, 04, 05, 06, 07, 08

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
% mRNAs14s{1} = data{7}.Data_sort(:,2); mu14s{1} = data{7}.mu;
% mRNAs14s{2} = data{6}.Data_sort(:,2); mu14s{2} = data{6}.mu;
% mRNAs14s{9} = data{3}.Data_sort(:,2); mu14s{9} = data{3}.mu;
% mRNAs14s{10} = data{4}.Data_sort(:,2); mu14s{10} = data{4}.mu;
mRNAs14s{11} = data{5}.Data_sort(:,2); mu14s{11} = data{5}.mu;
mRNAs14s{12} = data{8}.Data_sort(:,2); mu14s{12} = data{8}.mu;
x14s{11} = data{5}.Data_sort(:,1)/1000;
x14s{12} = data{8}.Data_sort(:,1)/1000;


slidedate = '2011-02-17/'; fname = 'MP05_22C_sna_y';  ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14s{1} = data{9}.Data_sort(:,2); mu14s{1} = data{9}.mu;
mRNAs14s{2} = data{12}.Data_sort(:,2); mu14s{2} = data{12}.mu;
mRNAs14s{3} = data{1}.Data_sort(:,2); mu14s{3} = data{1}.mu;
mRNAs14s{4} = data{2}.Data_sort(:,2); mu14s{4} = data{2}.mu;
mRNAs14s{5} = data{3}.Data_sort(:,2); mu14s{5} = data{3}.mu;
mRNAs14s{6} = data{5}.Data_sort(:,2); mu14s{6} = data{5}.mu;
mRNAs14s{13} = data{7}.Data_sort(:,2); mu14s{13} = data{7}.mu;
mRNAs14s{14} = data{8}.Data_sort(:,2); mu14s{14} = data{8}.mu;
x14s{1} = data{9}.Data_sort(:,1)/1000;
x14s{2} = data{12}.Data_sort(:,1)/1000;
x14s{3} = data{1}.Data_sort(:,1)/1000;
x14s{4} = data{2}.Data_sort(:,1)/1000;
x14s{5} = data{3}.Data_sort(:,1)/1000;
x14s{6} = data{5}.Data_sort(:,1)/1000;
x14s{13} = data{7}.Data_sort(:,1)/1000;
x14s{14} = data{8}.Data_sort(:,1)/1000;

slidedate = '2011-02-17/';   fname ='MP10_22C_sna_y_d'; ver = '_vN';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 

mRNAs14s{19} = data{4}.Data_sort(:,2); mu14s{19} = data{4}.mu(:,1);
mRNAs14s{20} = data{5}.Data_sort(:,2); mu14s{20} = data{5}.mu(:,1);
mRNAs14s{9} = data{2}.Data_sort(:,2); mu14s{9} = data{2}.mu(:,1);
mRNAs14s{10} = data{3}.Data_sort(:,2); mu14s{10} = data{3}.mu(:,1);
mRNAs14s{15} = data{7}.Data_sort(:,2); mu14s{15} = data{7}.mu(:,1);
mRNAs14s{16} = data{8}.Data_sort(:,2); mu14s{16} = data{8}.mu(:,1);
mRNAs14s{18} = data{1}.Data_sort(:,2); mu14s{18} = data{1}.mu(:,1);
x14s{19} = data{4}.Data_sort(:,1)/1000;
x14s{20} = data{5}.Data_sort(:,1)/1000;
x14s{9} = data{2}.Data_sort(:,1)/1000;
x14s{10} = data{3}.Data_sort(:,1)/1000;
x14s{15} = data{7}.Data_sort(:,1)/1000;
x14s{16} = data{8}.Data_sort(:,1)/1000;
x14s{18} = data{1}.Data_sort(:,1)/1000;

slidedate = '2011-02-17/'; fname = 'MP05_22C_sna_y_c';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14s{21} = data{1}.Data_sort(:,2); mu14s{21} = data{1}.mu(:,1);
mRNAs14s{22} = data{2}.Data_sort(:,2); mu14s{22} = data{2}.mu(:,1);
mRNAs14s{17} = data{4}.Data_sort(:,2); mu14s{17} = data{4}.mu(:,1);
x14s{21} = data{1}.Data_sort(:,1)/1000;
x14s{22} = data{2}.Data_sort(:,1)/1000;
x14s{17} = data{4}.Data_sort(:,1)/1000;

slidedate = '2011-12/'; fname = 's142_sna';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
 mRNAs14s{7} = data{14}.Data_sort(:,2); mu14s{7} = data{14}.mu;
 mRNAs14s{8} = data{16}.Data_sort(:,2); mu14s{8} = data{16}.mu;
x14s{7} = data{14}.Data_sort(:,1)/1000;
x14s{8} = data{16}.Data_sort(:,1)/1000;

end

% Combine data into cells sorted by age
mRNAsN{1} = 0; muN{1} = 0;  xN{1} = 0; 

mRNAs = [mRNAs12telo,mRNAs13,mRNAs13meta,mRNAs13telo,mRNAs14e,mRNAs14s,mRNAsN];
mu = [mu12telo,mu13,mu13meta,mu13telo,mu14e,mu14s,muN];

% keeping age grouping in separate cells
allx = {x12telo(2:end),x13,x13meta,x13telo,x14e,x14s,xN};
all_mRNA = {mRNAs12telo(2:end),mRNAs13,mRNAs13meta,mRNAs13telo,mRNAs14e,mRNAs14s,mRNAsN};
all_mu = {mu12telo(2:end),mu13,mu13meta,mu13telo,mu14e,mu14s,muN};



%% scatterhist

% get dimensions and initialize variables
G = length(all_mRNA);
Es = length(mRNAs);
mRNA = cell(Es,1); bkd=mRNA; 
CMap = zeros(Es,3);
% cmap = hsv(G); 
   
   cmap = [0,0,1;
        0,.75,1;
       .7,.4,1;
      .6,.2,.6;
       1,.2,.65;
      1,0,0;
       1,0,0];

Ls = [length(mRNAs12telo),length(mRNAs13),length(mRNAs13meta),length(mRNAs13telo),length(mRNAs14e),length(mRNAs14s),length(mRNAsN)];
g = [Ls(1),sum(Ls(1:2)),sum(Ls(1:3)),sum(Ls(1:4)),sum(Ls(1:5)),sum(Ls(1:6)),sum(Ls(1:7))];

group = 1;  % embryo group
gi = 0; % index of embryos within group

meds_on = zeros(1,Es); meds_off = zeros(1,Es); meds_offB = zeros(1,Es); 
Ncells = 0;
Nmolecules = 0;
Ncells_group = zeros(1,length(g));
Nmols_group = zeros(1,length(g));

dyn_fig = figure(2); clf; colordef white; set(gcf,'color','w'); hold on;
for e=2:Es
        
    if e > g(group)
        group = group + 1;
        gi = 0; 
    end
   % CMap(e,:) = [group/(length(g)-1),0,1-group/(length(g)-1)] ;
   CMap(e,:) = cmap(group,:); 
    
%         temp = sort(mu{e});
%         sna = mRNAs{e} - temp(2); %  min(mu{e});
%         pk = max(mu{e}-temp(2)); % min(mu{e}));
        sna = mRNAs{e} -  min(mu{e}(:,1));
        pk = max(mu{e}(:,1)- min(mu{e}(:,1)));
        mRNA{e} =  sna( sna>1/2*pk); 
        bkd{e} = sna( sna<1/2*pk);
        figure(2);  
        L = length(mRNA{e});
        plot(.7*rand(1,L)+e*ones(1,L),mRNA{e},'.','color',CMap(e,:),'MarkerSize',5);
        meds_on(e) = median(mRNA{e});
        
       Ncells = Ncells + length(sna);
       Nmolecules = Nmolecules + sum(sna); 
       Ncells_group(group) = Ncells_group(group)+ length(sna);
       Nmols_group(group) = Nmols_group(group)+ sum(sna);
       
        L2 = length(bkd{e});
        clr = CMap(e,:)*.15+[.85,.85,.85];
        while mean(clr)>=.9
            clr = clr -.05;
        end
        plot(.7*rand(1,L2)+e*ones(1,L2),bkd{e},'.','color',clr,'MarkerSize',5);
        % meds_off(e) = mean(bkd{e}(bkd{e}>0));
         meds_off(e) = sum(bkd{e}>20)/(L2+L);
          meds_offB(e) = sum(bkd{e}>20);  
end 
          
 ylim([0,300]); set(gca,'FontSize',14); xlim([0,75]);
ylabel('mRNA counts per cell');
xlabel(' "Time" (embryo number ascending age)');    
title([num2str(Nmolecules,'%10.0f') ' mRNAs from ',num2str(Ncells), ' cells']); 
%   %%  
    x = 2:Es;
%     
%     figure(2); hold on;
%     plot(x+.35,meds_on(2:end),'ko'); hold on;
%     plot(x+.35,smooth(meds_on(2:end)));
%     
%     figure(2);  hold on;
%     bkd_curve = smooth(meds_off(2:end));
%     plot(x+.35,250*bkd_curve); 
%     plot(x+.35,250*meds_off(2:end),'o');
    
%% neurogenic nuclei turn off (supp fig)
    figure(3); clf; colordef white; set(gcf,'color','w');
    bkd_curve = smooth(meds_off(2:end),.25);
    plot(x+.35,bkd_curve,'k-','linewidth',2); hold on;
    for e=2:Es
    plot(x(e-1)+.35,meds_off(e),'.','Color',CMap(e,:),'MarkerSize',25);
    end
    
%     bkd_curve = smooth(meds_offB(2:end),.25);
%     plot(x+.35,bkd_curve,'b-','linewidth',2); hold on;
%     plot(x+.35,meds_offB(2:end),'b+');
    
    
    ylim([0,.5]); set(gca,'FontSize',14);
    ylabel('fraction of partially active cells');
    
xlabel(' "Time" (embryo number, ascending age)');    
%% compute stats for dynamics calculations (synthesis and half-life)



mRNAon = cell(6,1);

ave_cnt = zeros(G,1); 
med_cnt = zeros(G,1);
upp_cnt = zeros(G,1);
low_cnt = zeros(G,1);
std_cnt = zeros(G,1); 
for k=1:G % k=3
    E = length(all_mRNA{k}) ;% length(all_mu{k})
    if k==3; E = 4; end
    for e=1:E 
        temp = sort(all_mu{k}{e});
%          sna = all_mRNA{k}{e} - temp(2);
%          pk = max(all_mu{k}{e}-temp(2) );
%         
            sna = all_mRNA{k}{e} - min(all_mu{k}{e}(:,1));
            pk = max(all_mu{k}{e}(:,1)-min(all_mu{k}{e}(:,1)));
            mRNAon{k}{e} =  sna(sna>1/2*pk); 
    end  
     med_cnt(k) = median(cellfun(@nanmean,mRNAon{k}));
     ave_cnt(k) = mean(cellfun(@nanmean,mRNAon{k}));
     
     sort_mean_cnt = sort(cellfun(@nanmean,mRNAon{k}));
     upp_cnt(k) = sort_mean_cnt(max([round(3/4*length(sort_mean_cnt)),1]));
     low_cnt(k) = sort_mean_cnt(max([1,round(1/4*length(sort_mean_cnt))]));
     std_cnt(k) = std(sort_mean_cnt);
end    

%%  half-life

% .5 = exp(-7.5/tau)
% tau = -7.5/log(.5) = 10.82
% exp(-4/10.82) 

prev_tau = 10.82;
T = 5; % estimate of mitosis length

m2 = mean(cellfun(@nanmean,mRNAon{4}));
m1 = mean(cellfun(@nanmean,mRNAon{3}));
m0 = mean(cellfun(@nanmean,mRNAon{1}));

%my_tau = -T/log( 2*med_cnt(4)/med_cnt(3) );
my_tau = -T/log( 2*m2/m1 )
max_tau =  max([-T/log( 2*upp_cnt(4)/med_cnt(3) ),-T/log( 2*med_cnt(4)/low_cnt(3))]);
min_tau = min([-T/log( 2*low_cnt(4)/med_cnt(3) ),-T/log( 2*med_cnt(4)/upp_cnt(3))]);

ubnd_tau = min([-T/log( 2*upp_cnt(4)/upp_cnt(3) ),-T/log( 2*upp_cnt(4)/upp_cnt(3))]);
lbnd_tau = min([-T/log( 2*low_cnt(4)/low_cnt(3) ),-T/log( 2*low_cnt(4)/low_cnt(3))]);

% tau = T/log(2m1/m2)
% 
% Dtau = sqrt( (dtau/dT*DT)^2 + dtau/dm1*Dm1)^2 + (dtau/dm2*Dm2)^2 )

% m1 = med_cnt(4); m2 = med_cnt(3); 
 
%  DT = .25; % 15 seconds uncertainty 
% Dtau = sqrt( (-1/log( 2*m1/m0)*DT)^2+(T/(log( 2*m1/m0)^2*m1)*std_cnt(4))^2 + (T/(log( 2*m1/m0)^2*m0)*std_cnt(3))^2)
% however error in tau is not semetric, so  

t = linspace(0,10,25);
degred = med_cnt(3)*exp(-t/my_tau); 
low_d = med_cnt(3)*exp(-t/max_tau);
high_d = med_cnt(3)*exp(-t/min_tau);

high_d = upp_cnt(3)*exp(-t/ubnd_tau);
low_d = low_cnt(3)*exp(-t/lbnd_tau); 
figure(4); clf; plot(t,degred,'k-','linewidth',3); hold on;
plot(t,low_d,'k--');

plot(0,med_cnt(3),'.','Color',cmap(3,:),'MarkerSize',30);
plot(T,2*med_cnt(4),'.','Color',cmap(4,:),'MarkerSize',30);
% plot(8,.5*med_cnt(3),'.','Color',[.6,.5,.1],'MarkerSize',30);
plot(-10,-10, 'k*','MarkerSize',10);
plot(0*ones(1,length(mRNAon{3})),cellfun(@mean,mRNAon{3}),'*','Color',cmap(3,:),'MarkerSize',10);
plot([T*ones(1,length(mRNAon{4}))],[2*cellfun(@mean,mRNAon{4})],'*','MarkerSize',10,'Color',cmap(4,:));
 plot(t,high_d,'k--');
 
xlim([-.5,7]);
legend(['measured snail halflife =',num2str(my_tau*log(2),3),' min ',...
    '(',num2str(min_tau*log(2),2),',',num2str(max_tau*log(2),2),')'],'error bounds',...
    'prometaphase','late telophase',...
    'mean counts from indiv. embryos','Location','Best');
set(gcf,'color','w');
ylabel('mRNA count','FontSize',15); xlabel('time (min)','FontSize',15);
title('mRNA degredation during cell division');
ylim([30,250]);



%% mRNA synthesis  with DNA replication

Tdiv = 14; % total length of division cycle 13 (McCleland et al 2008).  
Tstart = 2; % time for the first transcripts to start appearing
Trep = Tdiv/2; % time prior to locus replication
T2x = Trep - Tstart; %  time spent making mRNA prior to replication
T4x = Tdiv-Trep; % time after locus replication
t = linspace(0,Tdiv,100);
lifetime =  my_tau; % 1E10;% 
Tdiv_err = .2; % fractional uncertainty in timing or replication

% m(t) = s*tau*(1-exp(-T/tau))+m0*exp(-T/tau)

% m(t) = 2*s*tau*(1-exp(-t/tau))+m0*exp(-t/tau) 
%       + 2*s*tau*(1-exp(-(t)/tau))*H(t-T2)


m0 = ave_cnt(1); m1 =ave_cnt(3);  
synthesis_rate = (m1 - m0*exp(-(Tdiv)/lifetime))/ ( (2*lifetime*(1-exp(-(T2x+T4x)/lifetime)))+2*lifetime*(1-exp(-T4x/lifetime)) )  ;   
min_syn  = (low_cnt(3) - upp_cnt(1)*exp(-(Tdiv)/lifetime))/ ( (2*lifetime*(1-exp(-(T2x+T4x)/lifetime)))+2*lifetime*(1-exp(-T4x/lifetime)) )  ;   
max_syn = (upp_cnt(3) - low_cnt(1)*exp(-(Tdiv)/lifetime))/ ( (2*lifetime*(1-exp(-(T2x+T4x)/lifetime)))+2*lifetime*(1-exp(-T4x/lifetime)) )  ;   

tg1 = [zeros(1,ceil(100*Tstart/Tdiv)),linspace(0,T2x+T4x,100*(Tdiv-Tstart)/Tdiv)];
tg2 = [zeros(1,ceil(100*Trep/Tdiv)),linspace(0,T4x,100*(Tdiv-Trep)/Tdiv)];
prod = m0*exp(-t/lifetime) + 2*synthesis_rate*lifetime*(1-exp(-tg1/lifetime))  + 2*synthesis_rate*lifetime*(1-exp(-tg2/lifetime));

tgL = [zeros(1,ceil(100*Trep/Tdiv*(1+Tdiv_err))),linspace(0,T4x,100*(Tdiv-Trep)/Tdiv*(1-Tdiv_err))];
tgU = [zeros(1,ceil(100*Trep/Tdiv*(1-Tdiv_err))),linspace(0,T4x,100*(Tdiv-Trep)/Tdiv*(1+Tdiv_err))];
upp_prod =low_cnt(1)*exp(-t/lifetime) + 2*min_syn*lifetime*(1-exp(-tg1/lifetime))  + 2*min_syn*lifetime*(1-exp(-tgL/lifetime));
low_prod = upp_cnt(1)*exp(-t/lifetime) + 2*max_syn*lifetime*(1-exp(-tg1/lifetime))  + 2*max_syn*lifetime*(1-exp(-tgU/lifetime));


figure(5); clf; plot(t,prod,'k--','linewidth',3); hold on;
plot(t,upp_prod,'k--');
plot(0,m0,'.','color',cmap(1,:),'MarkerSize',30);
% plot(Tdiv/2,med_cnt(2),'.','Color',[2/8,0,1-2/8],'MarkerSize',30);
plot(Tdiv,m1,'.','Color',[4/8,0,1-4/8],'MarkerSize',30);
plot(-10,-10, 'k*','MarkerSize',10);

plot(0*ones(1,length(mRNAon{1})),cellfun(@mean,mRNAon{1}),'b*','MarkerSize',10);
plot(Tdiv*ones(1,length(mRNAon{3})),cellfun(@mean,mRNAon{3}),'*','Color',cmap(3,:),'MarkerSize',10);
plot(Tstart+1+linspace(0,Tdiv-Tstart-2,length(mRNAon{2})),sort(cellfun(@mean,mRNAon{2})),'*','Color',cmap(2,:),'MarkerSize',10);
 plot(t,low_prod,'k--');
set(gcf,'color','w'); %set(gca,'FontSize',15);
ylabel('mRNA count','FontSize',15); xlabel('time (min)','FontSize',15);

xlim([-.5,15]); ylim([50,250]);
mps = 60/synthesis_rate;
% ['tx rate =', num2str(synthesis_rate,2),'mRNA/min (', num2str(min_syn,2),',',num2str(max_syn,2), ')']

upp_mps = (1+sqrt( ((synthesis_rate - min_syn)/synthesis_rate).^2 + Tdiv_err.^2))*(mps);
low_mps =  (1-sqrt( ((synthesis_rate - max_syn)/synthesis_rate).^2 + Tdiv_err.^2))*(mps);
 
legend(['tx rate =', num2str(mps,2),'s/mRNA (', num2str(low_mps,2),',',num2str(upp_mps,2), ')'],...
    'uncertainty bound','ave # mRNA at telophase of cc12',...
    'ave # mRNA at prometaphase cc13',...
     'mean counts from indiv. embryos','Location','Best')
 title('mRNA synthesis rate');

 
 footprint = 50; elong_rate = 1000/60; 
 
 S_theory_max =   footprint/elong_rate ;
 
 %% predicted synthesis for cc14 without negative autoregulation
tg1 = [zeros(1,ceil(100*Tstart/Tdiv)),linspace(0,T2x+T4x,100*(Tdiv-Tstart)/Tdiv)];
tg2 = [zeros(1,ceil(100*Trep/Tdiv)),linspace(0,T4x,100*(Tdiv-Trep)/Tdiv)];
prod = m0*exp(-t/lifetime) + 4*synthesis_rate*lifetime*(1-exp(-tg1/lifetime))  + 4*synthesis_rate*lifetime*(1-exp(-tg2/lifetime));

 figure(14); clf; plot(t,prod); title('predicted cc13 4x synthesis'); 
 
 
%% predicted synthesis for cc14 without negative autoregulation
Trep = 10;
Tdiv = 60;
T2x = Trep - Tstart; %  time spent making mRNA prior to replication
T4x = Tdiv-Trep; % time after locus replication
lifetime =  my_tau; % 1E10; %  
t = linspace(0,Tdiv,100);
tg1 = [zeros(1,ceil(100*Tstart/Tdiv)),linspace(0,T2x+T4x,100*(Tdiv-Tstart)/Tdiv)];
tg2 = [zeros(1,ceil(100*Trep/Tdiv)),linspace(0,T4x,100*(Tdiv-Trep)/Tdiv)];
prod = m0*exp(-t/lifetime) + 2*synthesis_rate*lifetime*(1-exp(-tg1/lifetime))  + 2*synthesis_rate*lifetime*(1-exp(-tg2/lifetime));

 figure(15); clf; plot(t,prod); title('predicted cc14 synthesis'); 
  
%%
fout = 'C:\Users/Alistair/My Documents/Projects/Snail Patterning/Results/';
if savedata == 1
saveas(dyn_fig,[fout,'snail_dyn.eps'],'eps');
end
%% plot representative nuclei from each staging point
if loadnucdata == 1
    fout = 'C:\Users/Alistair/My Documents/Projects/Snail Patterning/Results/';
    Nfig = figure(3); clf; 
    e = 2; Nchn = Nuc{e}(:,:,3); 
    imagesc(Nchn); colormap copper; 
    caxis([mean(Nchn(:))/1.1,3*mean(Nchn(:))]); 
    set(gcf,'color','k');
    if savedata==1
    saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);
    end

    Nfig = figure(3); clf; 
    e = 10; Nchn = Nuc{e}(:,:,3); 
    imagesc(Nchn); colormap copper; 
    caxis([mean(Nchn(:))/2,2*mean(Nchn(:))]); 
    set(gcf,'color','k');
    if savedata==1
    saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);
    end

    Nfig = figure(3); clf; 
    e = 13; Nchn = Nuc{e}(:,:,3); 
    imagesc(Nchn); colormap copper; 
    caxis([mean(Nchn(:))/2,3*mean(Nchn(:))]); 
    set(gcf,'color','k');
    if savedata==1
    saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);
    end

    Nfig = figure(3); clf;
    e = 15; Nchn = Nuc{e}(:,:,3); 
    imagesc(Nchn); colormap copper; 
    caxis([mean(Nchn(:))/2,4*mean(Nchn(:))]); 
    set(gcf,'color','k');
    if savedata==1
    saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);
    end

    Nfig = figure(3); clf;
    e = 21; Nchn = Nuc{e}(:,:,3); 
    imagesc(Nchn); colormap copper; 
    caxis([mean(Nchn(:))/4,4*mean(Nchn(:))]);
    set(gcf,'color','k');
    if savedata==1
    saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);
    end

    Nfig = figure(3); clf;
    e = 23; Nchn = Nuc{e}(:,:,3); 
    imagesc(Nchn); colormap copper; 
    caxis([mean(Nchn(:))/4,4*mean(Nchn(:))]); 
    set(gcf,'color','k');
    if savedata==1
    saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);
    end

    Nfig = figure(3); clf;
    e = 28; Nchn = Nuc{e}(:,:,3);  % 28
    imagesc(Nchn); colormap copper; 
    caxis([mean(Nchn(:))/4,4*mean(Nchn(:))]); 
    set(gcf,'color','k');
    if savedata==1
    saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);
    end


    Nfig = figure(3); clf;
    e = 56; Nchn = Nuc{e}(:,:,3);  % 28
    imagesc(Nchn); colormap copper; 
    caxis([mean(Nchn(:))/1.5,2*mean(Nchn(:))]); 
    set(gcf,'color','k');
    if savedata==1
    saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);
    end

end




%% boxplot mRNA, CoV and Fano

labs = {'cc12telo','cc13int','cc13meta','cc13telo','cc14early','cc14mid',''};
aveon = NaN*ones(G,50);
cov = NaN*ones(G,50);
fano = NaN*ones(G,50);

F = 15; % fontsize;


for k=1:G % k=3
    for e=1:length(all_mRNA{k}) 
            sna = all_mRNA{k}{e} - min(all_mu{k}{e}(:,1));
            pk = max(all_mu{k}{e}(:,1)-min(all_mu{k}{e}(:,1)));
            tmp = sna(sna>1/2*pk);
            aveon(k,e) = nanmean(tmp);
            cov(k,e) = nanstd(tmp)/ nanmean(tmp);
            fano(k,e) = nanstd(tmp)^2/ nanmean(tmp);
    end  
end    

figure(10); clf; boxplot(aveon','Labels',labs,'colors',cmap, 'plotstyle','compact','whisker',10);
ylabel('mRNA count','FontSize',F);
set(gcf,'color','w'); set(gca,'color','w','FontSize',F);
hand = findobj(gca,'linewidth',4); set(hand,'linewidth',20);
whisk = findobj(gca,'-not','color','k','-and','linewidth',.5);
set(whisk,'LineStyle','--','linewidth',2);

set(gca,'color','w','FontSize',F); set(gcf,'color','w');
figure(11); clf; boxplot(cov','Labels',labs,'colors',cmap, 'plotstyle','compact','whisker',10); ylabel('CoV','FontSize',F);
set(gcf,'color','w'); set(gca,'color','w','FontSize',F);
hand = findobj(gca,'linewidth',4); set(hand,'linewidth',20);
whisk = findobj(gca,'-not','color','k','-and','linewidth',.5);
set(whisk,'LineStyle','--','linewidth',2);

figure(12); clf; boxplot(fano','Labels',labs,'colors',cmap, 'plotstyle','compact','whisker',10); ylabel('Fano factor','FontSize',F);
set(gcf,'color','w'); set(gca,'color','w','FontSize',F);
hand = findobj(gca,'linewidth',4); set(hand,'linewidth',20);
whisk = findobj(gca,'-not','color','k','-and','linewidth',.5);
set(whisk,'LineStyle','--','linewidth',2);

%% spatial curves

addpath('C:\Users\Alistair\Documents\Projects\GenCode');

    fit_sigs = 1;     xmin = 0; xmax = 1; 
    hillcoefs = zeros(50,G);
    offsets  = zeros(50,G);
    sigmax = zeros(50,G);
    bkd  = zeros(50,G);
    fits = cell(50,G);
    celldat = cell(50,G);
  % align gradients 
  for g=1:G-1
      E = length(all_mRNA{g});
      for e=1:E
               if fit_sigs == 1 % Sigmoidal fitting routine
                  xdat = allx{g}{e}; % renaming data, convert to microns
                  sna = all_mRNA{g}{e} - min(all_mu{g}{e}(:,1)); % subtract background 
                  grad = sna; L = length(grad); st = floor(L/5);
                 
                  % truncate data from right to get only sigmoids and not full steps  
                 grad = grad(st:end); 
                 gx = xdat(st:end);
                 n = 3; % initial guesses for sigmoid fits
                 theta = 2*mean(gx);   % threshold
                 A = max(grad); % fxn max
                 b=min(grad);   % offset / background
                 [p,fit] = fxn_fit_sigmoid(gx',grad',[n,theta,A,b],'r');
                 hillcoefs(e,g) = p(1);
                 offsets(e,g) = abs(p(2));
                 sigmax(e,g) = p(3); 
                 bkd(e,g) = p(4); 
                 fits{e,g} = [gx'; fit]; 
                 celldat{e,g} = [xdat';sna'];
                 
                 xmin = min(xmin,-offsets(e,g)); % for global plotting
                 xmax = max(xmax, max(xdat)); % for global plotting
               end         
      end
  end
  
%% Plot curves  
figure(21); clf; figure(22); clf;
pts = 1000; % number of points for interpolation
stp = 20; % interval at which to plot spread bars
xs = linspace(xmin,xmax-xmin,pts); % for interp. must be positive
xo = xs + xmin; % for plotting, 0 is located at snail 'boundary'

for g=1:G-1
    E = length(allx{g});
    c = zeros(E,pts); 
    for e = 1:E
        if isnan(offsets(e,g)) == 0 % this is NaN if embryo e is not in filtered datatset 
             x = celldat{e,g}(1,:)-offsets(e,g) - xmin;
             y = celldat{e,g}(2,:);
             [x,v] = unique(x); % required for interp to work
             y = y(v);
             c(e,:) = interp1(x',y',xs); 
%     figure(21); plot(x+xmin,y,'o','color',CMap(g,:),'MarkerSize',1); 
%     hold on; plot(xs+xmin,c(e,:),'.','color',CMap(g,:),'MarkerSize',1);
        end
    end
        c(c==0) = NaN;
        cave = nanmean(c);
        cstd = nanstd(c);
        figure(22);  plot(xo,cave,'color',cmap(g,:),'linewidth',3); hold on;
        plot(xo,cave-.5*cstd,'--','color',cmap(g,:)); hold on;
        plot(xo,cave+.5*cstd,'--','color',cmap(g,:)); hold on;
              
        X = [xo,xo]; Y = [cave-.5*cstd,cave+.5*cstd];
        Y(isnan(Y)) = 0; % can't have NaN's in plot 
        
        verts = [X',Y']; % 
        faces = [linspace(1,pts,pts),fliplr(linspace(pts+1,2*pts,pts))];      
        
        figure(21); 
        p = patch('Faces',faces,'Vertices',verts,'FaceColor',cmap(g,:),'EdgeColor','none');
        hold on;  alpha .7;
         plot(xo,cave,'color',cmap(g,:),'linewidth',2); hold on;
        
end
figure(22);
   xlim([-4E1,3E1]); xlabel('distance from boundary, (um)','FontSize',F);
   ylim([0,1.2*max(aveon(:))]);  ylabel('mRNA count','FontSize',F); 
   set(gca,'FontSize',F); % legend(groupnames); 
   set(gcf,'color','w');
   
   figure(21); 
   xlim([-4E1,3E1]); xlabel('distance from boundary, (um)','FontSize',F);
   ylim([0,1.2*max(aveon(:))]);  ylabel('mRNA count','FontSize',F); 
   set(gca,'FontSize',F); % legend(groupnames); 
      set(gcf,'color','w');
   
  addpath('../../GenCode/');

                