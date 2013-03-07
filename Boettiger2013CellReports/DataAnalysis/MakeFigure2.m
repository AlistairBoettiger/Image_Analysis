%%                          Fig_snail_vs_UASlacZ.m
% Alistair Boettiger                                   Date Begun: 05/30/12
% Zhuang Lab / Levine Lab                           Last Modified: 07/11/12
clear all;

% adapted from snail dynamics figure
% compare CoVs 

%% Publication figure as of 09/26/12

%%
% load all data sorted by age

    folder = 'C:\Users/Alistair/My Documents/Projects/mRNA_counting/Data/';
    rawfolder = 'D:/Data/';    
    MGafolder = 'C:\Users\Alistair\Documents\Projects\mRNA_counting\Data\2011-02-17\';
    load([MGafolder,'Mga_data']); % Mga mean and std expression data 
  snail_results ='C:\Users\Alistair\Documents\Projects\Snail Patterning\Results\';
  
   savedata = 0; % don't write data to disk (will auto-overwrite existing)
addpath('C:\Users/Alistair/My Documents/Projects/GenCode');
 
   
 %% Load mRNA data, sort by age  
 % data is hand-sorted into age categories.  Should do this automatically
 % during original image processing in the future.   
 clear mRNAs mu; 
 
% ----------------------------% exiting cc12 ---------------------------%
% 'wt_sna'  emb = 03, 04, 10
% 'MP06_cflip' emb = 08
% 'MP06Hz_b' emb = 02
slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{5} = data{3}.Data_sort(:,2);  mu12telo{5} = data{3}.mu;
mRNAs12telo{3} = data{4}.Data_sort(:,2);  mu12telo{3} = data{4}.mu;
% mRNAs12telo{1} = data{10}.Data_sort(:,2); mu12telo{1} = data{10}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{4} = data{8}.Data_sort(:,2); mu12telo{4} = data{8}.mu;
%mRNAs{4} = data{8}.Data_sort(:,3); mu{4} = data{8}.mu(:,2); vout = '_o2'

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{2} = data{2}.Data_sort(:,2); mu12telo{2} = data{2}.mu;

slidedate = '2011-12/'; fname = 's142_sna';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{6} = data{1}.Data_sort(:,2); mu12telo{6} = data{1}.mu;
mRNAs12telo{7} = data{2}.Data_sort(:,2); mu12telo{7} = data{2}.mu;

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

slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{2} = data{9}.Data_sort(:,2); mu13{2} = data{9}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{3} = data{1}.Data_sort(:,2); mu13{3} = data{1}.mu;
mRNAs13{4} = data{2}.Data_sort(:,2); mu13{4} = data{2}.mu;
mRNAs13{10} = data{9}.Data_sort(:,2); mu13{10} = data{9}.mu;
mRNAs13{8} = data{10}.Data_sort(:,2); mu13{8} = data{10}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{6} = data{4}.Data_sort(:,2); mu13{6} = data{4}.mu;
mRNAs13{7} = data{5}.Data_sort(:,2); mu13{7} = data{5}.mu;
%========================================================================%



% --------------------- % metaphase cc13 ----------------------------------
% 'wt_sna'  emb = 05, 06 
% 'MP06_cflip' emb = 03, 04, 06
% 'MP06Hz_b' emb = 01, 09, 10

slidedate = '2011-12/'; fname = 's142_sna';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{1} = data{15}.Data_sort(:,2); mu13meta{1} = data{15}.mu;

slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{4} = data{5}.Data_sort(:,2); mu13meta{4} = data{5}.mu;
mRNAs13meta{3} = data{6}.Data_sort(:,2); mu13meta{3} = data{6}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{2} = data{3}.Data_sort(:,2); mu13meta{2} = data{3}.mu;
mRNAs13meta{9} = data{4}.Data_sort(:,2); mu13meta{9} = data{4}.mu;
% mRNAs13meta{6} = data{6}.Data_sort(:,2); mu13meta{5} = data{6}.mu;

slidedate = '2011-12/';   fname ='s140_sna'; ver = '_v2';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data');
mRNAs13meta{6} = data{7}.Data_sort(:,2); mu13meta{6} = data{7}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{7} = data{1}.Data_sort(:,2); mu13meta{7} = data{1}.mu;
mRNAs13meta{8} = data{9}.Data_sort(:,2); mu13meta{8} = data{9}.mu;
mRNAs13meta{5} = data{10}.Data_sort(:,2); mu13meta{5} = data{10}.mu;
%------------------------------------------------------------------%

% =============== % telophase into cc14 ===============================
% 'MP06_cflip' emb = 05 
% 'MP06Hz_b' emb = 03, 07
% 'MP06Hz' emb = 10
slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{3} = data{5}.Data_sort(:,2); mu13telo{3} = data{5}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{4} = data{3}.Data_sort(:,2); mu13telo{4} = data{3}.mu;
mRNAs13telo{1} = data{7}.Data_sort(:,2); mu13telo{1} = data{7}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{2} = data{10}.Data_sort(:,2);  mu13telo{2} = data{10}.mu;
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

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14e{5} = data{6}.Data_sort(:,2);   mu14e{5} = data{6}.mu;
mRNAs14e{14} = data{8}.Data_sort(:,2);  mu14e{14} = data{8}.mu;
mRNAs14e{9} = data{11}.Data_sort(:,2);  mu14e{9} = data{11}.mu;


slidedate = '2011-12/';   fname ='s140_sna'; ver = '_v2';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data');
mRNAs14e{12} = data{1}.Data_sort(:,2); mu14e{12} = data{1}.mu(:,1);
mRNAs14e{3} = data{3}.Data_sort(:,2); mu14e{3} = data{3}.mu(:,1);
mRNAs14e{15} = data{4}.Data_sort(:,2); mu14e{15} = data{4}.mu(:,1);
mRNAs14e{6} = data{5}.Data_sort(:,2); mu14e{6} = data{5}.mu(:,1);
mRNAs14e{13} = data{8}.Data_sort(:,2); mu14e{13} = data{8}.mu(:,1);
mRNAs14e{16} = data{9}.Data_sort(:,2); mu14e{16} = data{9}.mu(:,1);

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


slidedate = '2011-02-17/';   fname ='MP10_22C_sna_y_d'; ver = '_vN';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 

mRNAs14s{19} = data{4}.Data_sort(:,2); mu14s{19} = data{4}.mu(:,1);
mRNAs14s{20} = data{5}.Data_sort(:,2); mu14s{20} = data{5}.mu(:,1);
mRNAs14s{9} = data{2}.Data_sort(:,2); mu14s{9} = data{2}.mu(:,1);
mRNAs14s{10} = data{3}.Data_sort(:,2); mu14s{10} = data{3}.mu(:,1);
mRNAs14s{15} = data{7}.Data_sort(:,2); mu14s{15} = data{7}.mu(:,1);
mRNAs14s{16} = data{8}.Data_sort(:,2); mu14s{16} = data{8}.mu(:,1);
mRNAs14s{18} = data{1}.Data_sort(:,2); mu14s{18} = data{1}.mu(:,1);

slidedate = '2011-02-17/'; fname = 'MP05_22C_sna_y_c';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14s{21} = data{1}.Data_sort(:,2); mu14s{21} = data{1}.mu(:,1);
mRNAs14s{22} = data{2}.Data_sort(:,2); mu14s{22} = data{2}.mu(:,1);
mRNAs14s{17} = data{4}.Data_sort(:,2); mu14s{17} = data{4}.mu(:,1);

slidedate = '2011-12/'; fname = 's142_sna';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
 mRNAs14s{7} = data{14}.Data_sort(:,2); mu14s{7} = data{14}.mu;
 mRNAs14s{8} = data{16}.Data_sort(:,2); mu14s{8} = data{16}.mu;

% combine data in cell arrays by age
mRNAsN{1} = 0; muN{1} = 0; 
mRNAs = [mRNAs12telo,mRNAs13,mRNAs13meta,mRNAs13telo,mRNAs14e,mRNAs14s,mRNAsN];
mu = [mu12telo,mu13,mu13meta,mu13telo,mu14e,mu14s,muN];

all_mRNA = {mRNAs12telo(2:end),mRNAs13,mRNAs13meta,mRNAs13telo,mRNAs14e,mRNAs14s,mRNAsN};
all_mu = {mu12telo(2:end),mu13,mu13meta,mu13telo,mu14e,mu14s,muN};
G = length(all_mRNA);

% save([snail_results,'mRNA_cnt_data.mat']);

%% boxplot dynamics: mRNA, CoV, Fano
labs = {'cc12telo','cc13int','cc13meta','cc13telo','cc14early','cc14mid',''};
aveon = NaN*ones(G,50);
cov = NaN*ones(G,50);
fano = NaN*ones(G,50);

F = 15; % fontsize;

  cmap = [0,0,1;
        0,.75,1;
       .7,.4,1;
      .6,.2,.6;
       1,.2,.65;
      1,0,0;
       1,0,0];

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

LacZ_cnt = [MGa1x_ave(:,1);MGa2x_ave(:,1)];
LacZ_std = [MGa1x_std(:,1);MGa1x_std(:,1)];
LacZ_fano = LacZ_std.^2./LacZ_cnt;
LacZ_cov = LacZ_std./LacZ_cnt;

N_lacZ = length(LacZ_cnt);
N_sna = length(all_mRNA{6});
N_sna_Tcells = 6963;
N_lacZ_Tcells = 4845;

comp_fano = [ fano(6,:)', [LacZ_fano; NaN*ones(50-N_lacZ,1)]];
comp_cnt = [ aveon(6,:)',[LacZ_cnt; NaN*ones(50-N_lacZ,1)]];
comp_cov = [ cov(6,:)', [LacZ_cov; NaN*ones(50-N_lacZ,1)]];

labs = {['D. mel. sna N=',num2str(N_sna),' embryos, ',num2str(N_sna_Tcells),' cells'] ;
       ['D. mel. UAS-LacZ N=',num2str(N_lacZ),' embryos, ',num2str(N_lacZ_Tcells),' cells']} ;
figure(10); clf; boxplot(comp_fano,'Labels',labs,'colors',cmap([1,5],:), 'plotstyle','compact','whisker',10); 
ylabel('Fano factor','FontSize',F); ylim([0,30]);
set(gcf,'color','w'); set(gca,'color','w','FontSize',F);
hand = findobj(gca,'linewidth',4); set(hand,'linewidth',20);
whisk = findobj(gca,'-not','color','k','-and','linewidth',.5);
set(whisk,'LineStyle','--','linewidth',2);

figure(11); clf; boxplot(comp_cnt,'Labels',labs,'colors',cmap([1,5],:), 'plotstyle','compact','whisker',10); 
ylabel('mRNA per cell','FontSize',F); ylim([0,250]);
set(gcf,'color','w'); set(gca,'color','w','FontSize',F);
hand = findobj(gca,'linewidth',4); set(hand,'linewidth',20);
whisk = findobj(gca,'-not','color','k','-and','linewidth',.5);
set(whisk,'LineStyle','--','linewidth',2);

figure(12); clf; boxplot(comp_cov,'Labels',labs,'colors',cmap([1,5],:), 'plotstyle','compact','whisker',10); 
ylabel('CoV','FontSize',F); ylim([0,1]);
set(gcf,'color','w'); set(gca,'color','w','FontSize',F);
hand = findobj(gca,'linewidth',4); set(hand,'linewidth',20);
whisk = findobj(gca,'-not','color','k','-and','linewidth',.5);
set(whisk,'LineStyle','--','linewidth',2);


%%  Compare CoV betweeen Species and Studies

med_CoV = nanmedian(cov'); 
sna_dyn_cnts = nanmedian(aveon');
sna_stds = nanstd(aveon');

    addpath('C:\Users\Alistair\Documents\Projects\GenCode\');
    cc14l_err = bserr(cov(6,:),'nanmedian');
    cc14t_err = bserr(cov(3,:),'nanmedian');
    Lacz1_err = bserr(eta1x,'nanmedian');
    Lacz2_err = bserr(eta2x,'nanmedian');
    
    scr_cov = 32.38/94.47;
    scr_err = std([26.49/99.56,29.69/104.4,35.89/80]);
  
     
   
    
  sna_cnts = [sna_dyn_cnts(1), sna_dyn_cnts(6), nanmedian(MGa1x_ave(:,1)), nanmedian(MGa2x_ave(:,1)), 94.47];
  sna_stds = [sna_stds(1), sna_stds(6), nanmedian(MGa1x_std(:,1)), nanmedian(MGa2x_std(:,1)), 32.38];
  sna_COV =  [med_CoV(1),
              med_CoV(6), 
              nanmedian( MGa1x_std(:,1)./MGa1x_ave(:,1) ),
              nanmedian( MGa2x_std(:,1)./MGa2x_ave(:,1) ),
              scr_cov]';  
  sna_err = [ cc14l_err; cc14t_err; Lacz1_err; Lacz2_err; scr_err];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~ Yeast measurements ~~~~~~~~~~~~~~~~~~~~~~~~~~%
% From Zenklusen et al 2008 NSMB
  MDN1 = []; f1 = 1000*[0, .02, .04, .09, .12, .17, .15, .11, .12, .07, .04, .03, .02, .01, 0 ,0,zeros(1,20)];
  KAP104 = []; f2 = 1000*[.02,.025, .09, .14, .23, .14, .1, .11, .1, .06, .002,.01, .005,0,0,0,zeros(1,20)]; 
  DOA1 = []; f3 = 1000*[.09,.19,.28,.2,.12,.09,.05,.01,.005,0,0,0,0,0,0,0,zeros(1,20)];
  POL1 = []; f4 = 1000*[.23,.23,.14,.1,.09,.04,.03,.04,.045,.03,.04,.015,.025,.01,0,0,zeros(1,20)];
  PDR5 = []; f5 = 1000*[.01 .01 .02 .035 .035 .05 .05 .065 .065 .05 .05 .06 .06 .05 .025 .04 .04 .02 .02 .025 .03 .02 .01 .02 .02 .02 .01 .001 .02 .02 0 .01 .01 .01];
  for nt = 1:34
      MDN1 = [MDN1, (nt-1)*ones(1,f1(nt))];
      KAP104 = [KAP104, (nt-1)*ones(1,f2(nt))];
      DOA1 = [DOA1, (nt-1)*ones(1,f3(nt))];
      POL1 = [POL1, (nt-1)*ones(1,f4(nt))];
      PDR5 = [PDR5, (nt-1)*ones(1,f5(nt))];
  end
  figure(1); clf; 
  subplot(5,1,1); hist(MDN1,[0:15]);
  subplot(5,1,2); hist(KAP104,[0:15]);
  subplot(5,1,3); hist(DOA1,[0:15]);
  subplot(5,1,4); hist(POL1,[0:15]);
  subplot(5,1,5); hist(PDR5,[0:40]);

   % last two data-points from To et al 2010 Science
 yeast_cnts = [mean(MDN1),mean(KAP104),mean(DOA1),mean(POL1),mean(PDR5),2.8,3.7]; % 0 data not available 
 yeast_stds = [std(MDN1),std(KAP104),std(DOA1),std(POL1),std(PDR5),2.8,1.15*3.7];
 yeast_COV = yeast_stds./yeast_cnts; 
%  
% yeast_COV = [std(MDN1)/mean(MDN1),
%    std(KAP104)/mean(KAP104),
%    std(DOA1)/mean(DOA1),
%    std(POL1)/mean(POL1),
%    std(PDR5)/mean(PDR5),
%    2.5, 
%    1.1];

% % Fano Factor
%     var(MDN1)/mean(MDN1)
%     var(KAP104)/mean(KAP104)
%     var(DOA1)/mean(DOA1)
%     var(POL1)/mean(POL1)
%     var(PDR5)/mean(PDR5)
  yeast_err = [.1 .1 .1 .2 .24 .2 .1]';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% ~~~~~~~~~~~~~~~ Worm CoV (Raj 2010) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% These are end-1 for different genotypes!  
% 'N2';'zu67';'zu135';'zu129';
worm_COV =  [.23 ;   .33;  .29;  .38];
  worm_err = [.03 .04 .07  .05]';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%~~~~~~~~~~~~~~~~~~~ COS Cell (Raj 2006) ~~~~~~~~~~~~~~~~~~~~~~~~~%
% (see excell document)
Mcell_cnts = [42.05, 34.83, 78.80, 41.71]; 
Mcell_stds = [44.86, 35.98, 145.35, 122.20];
Mcell_COV = Mcell_stds./Mcell_cnts; 
Mcell_err = [.18 .1 .25 .75]';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


% ~~~~~~ Mouse Intenstine CoV (Ishkovitz Nature Cell Bio 2011) ~~~~~~~~~~%
%   % Estimates by eye from supp mat
%   % Prom-1 Bmpr1a EphB2 CD44 Bmi1 Ascl2 Msi-1 Olfm4 EphB3, mTert, Sox9,
%   % sox4, Lgr5, mmp7, DCAMKL1
% mouse_COV = [.6 .7 .8 1 1.08 1.2 1.25 1.42 1.45 1.5 1.7 1.8 1.85 2.15 3.92]';
% mouse_err = [0  0   0  0  0   0    0    0   0    0    0  0   0    0    0]'; 
 % clear all;

load([snail_results,'Itzkovitz_mouse2.mat'])

j = 0; genes = cell(30*3,1); 
mouse_stem = cell(30*3,1); 
mouse_dif = cell(30*3,1); 

for i=1:30;
    for c=1:3
        j = j+1;
        s_gene = strfind(tag_full_cont{i}{c},'-');
        e_gene = strfind(tag_full_cont{i}{c},'_');
        genes{j} = tag_full_cont{i}{c}(s_gene+1:e_gene-1);
        mouse_stem{j} = profile_full_abs_cont{i}(:,c+1);
        mouse_dif{j} = profile_full_abs_cont{i}(:,c+1);
    end
end
ncells_tot = cellfun(@length,mouse_stem); % total cells

j=0;
n_genes = unique(genes);
mouse_stem_cnt = cell(1,length(n_genes));
mouse_dif_cnt = cell(1,length(n_genes));


for i=1:30
    for c=1:3
        j = j+1;
        geneid =  strmatch(genes{j},n_genes);
        filt =  abs(profile_full_abs_cont{i}(:,1))<5;
        mouse_stem_cnt{geneid} = [mouse_stem_cnt{geneid}; mouse_stem{j}(filt)];
        filt =  abs(profile_full_abs_cont{i}(:,1))>4;
        mouse_dif_cnt{geneid} = [mouse_dif_cnt{geneid}; mouse_dif{j}(filt)];
    end
end

ncells_filt = cellfun(@length,mouse_stem_cnt) % filtered 'comparable' cells

m_stem_means = cellfun(@mean,mouse_stem_cnt);
m_stem_stds = cellfun(@std,mouse_stem_cnt);
m_stem_COV = m_stem_stds./m_stem_means
m_dif_means = cellfun(@mean,mouse_dif_cnt);
m_dif_stds = cellfun(@std,mouse_dif_cnt);
m_dif_COV = m_dif_stds./m_dif_means

% m_stem_COV(m_stem_means<1) = NaN
% m_dif_COV(m_dif_means<1) = NaN

mouse_err = zeros(length(n_genes),1); 
for j=1:length(n_genes)
    mouse_err(j) = (bserr(mouse_stem_cnt{j},'std'))/ (m_stem_means(j));
end
mouse_COV = m_stem_COV';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% figure(11); clf; plot([sna_COV; worm_COV; yeast_COV; Mcell_COV],'k.');
% figure(11); clf; plot(log([sna_COV; worm_COV; yeast_COV; Mcell_COV]),'ko');
  
   labs = {'sna early';'sna steady state';'LacZ-1'; 'LacZ-2'; 'Scr (PS2)';...
        'Asc12'; 'Bmi1'; 'CD44'; 'DCAMKL1'; 'LGRX'; 'Msi1'; 'Olfm4'; 'P53';
        'MDN1';'KAP104';'DOA1';'POL1';'PDR5';'7xTetO';'1xTetO';...
        '1x Tet no dox';'1x Tet low dox';'7x Tet no dox';'7x Tet low dox'};
 
    empty_labs = {'  ';'  ';'  '; '  '; '  ';...
        '  '; '  '; '  '; '  '; '  '; '  '; '  '; '  ';
        '  ';'  ';'  ';'  ';'  ';'  ';'  ';...
        '  ';'  ';'  ';'  '};
 
    
 % Plotting    
    Cmap = [];
    for i=1:length(sna_COV)
        Cmap = [Cmap; cmap(1,:)];
    end
    for i=1:length(mouse_COV)
        Cmap = [Cmap; cmap(2,:)];
    end
    for i=1:length(yeast_COV)
        Cmap = [Cmap; cmap(4,:)];
    end
    for i=1:length(Mcell_COV)
        Cmap = [Cmap;cmap(6,:)];
    end  
    
    all_cnts = [sna_cnts,  m_stem_means, yeast_cnts, Mcell_cnts];
    all_cnt_std =[sna_stds, m_stem_stds,yeast_stds, Mcell_stds];
    all_COV = [ sna_COV'; mouse_COV; yeast_COV'; Mcell_COV' ];
    all_err = [ sna_err; mouse_err; yeast_err; Mcell_err  ];
    
    
  figure(1); clf; set(gcf,'color','w'); colordef white;
  subplot(4,1,1);
     barweb(all_COV,all_err,[],'','','','',Cmap);
     set(gca,'Ytick',[1,4],'FontSize',15);
     ylim([.9,4]);
   xticklabel_rotate(linspace(.61,1.38,length(labs)),70,empty_labs,'interpreter','none','FontSize',16)
subplot(4,1,2:4);
     barweb(all_COV,all_err,[],[],[],[],[],Cmap); 
     set(gca,'Ytick',[.1,.3,.5,.7,.9],'FontSize',15);
     ylim([.08,.9]);
     ylabel('Cov','FontSize',15);
   xticklabel_rotate(linspace(.61,1.38,length(labs)),70,labs,'interpreter','none','FontSize',6)
   
   
     figure(2); clf; colordef white;
     barweb(all_COV,all_err,[],'','','','',Cmap);
     set(gca,'Ytick',[1,2,3,4],'FontSize',15); 
     ylim([.08,4]); ylabel('CoV','FontSize',15);
      xticklabel_rotate(linspace(.61,1.38,length(labs)),70,labs,'interpreter','none','FontSize',10);
      set(gcf,'color','w');
    
    % mRNA Counts
      figure(5); clf; colordef white;
    % my errorbars for log plots:  
      v1 = log10(all_cnts); 
      x1 = linspace(.615,1.383,length(labs)); x2 = [];  xgap = x1(2)-x1(1);
      e1 = log10(all_cnts+all_cnt_std); e2 = [];
      for n=1:length(labs)
          x2 = [x2,x1(n),x1(n),NaN,x1(n)-xgap/2,x1(n)+xgap/2,NaN];
          e2 = [e2,v1(n),e1(n),NaN,e1(n),e1(n),NaN];
      end
     hold on; plot(x2,e2,'k-','linewidth',2);
     barweb(log10(all_cnts),zeros(size(all_cnt_std)),[],'','','','',Cmap, [], [], 1);
     set(gca,'FontSize',16); 
     ylabel('log_{10}(mRNA per cell)','FontSize',16);
      xticklabel_rotate(linspace(.61,1.38,length(labs)),70,labs,'interpreter','none','FontSize',10);
      set(gcf,'color','w');
 
     
      
      % Fano Factors 
      fano_err = zeros(size(all_cnts)); % log10((all_err'.*all_cnts).^2./all_cnts)
      figure(6); clf; colordef white;
          % my errorbars for log plots:  
      v1 = log10(all_cnt_std.^2./all_cnts); 
      x1 = linspace(.615,1.383,length(labs)); x2 = [];  xgap = x1(2)-x1(1);
      e1 = log10(all_cnt_std.^2./all_cnts+(all_err'.^2).*all_cnts); e2 = [];
      for n=1:length(labs)
          x2 = [x2,x1(n),x1(n),NaN,x1(n)-xgap/2,x1(n)+xgap/2,NaN];
          e2 = [e2,v1(n),e1(n),NaN,e1(n),e1(n),NaN];
      end
     hold on; plot(x2,e2,'k-','linewidth',2);  
     barweb(v1, zeros(size(all_cnts)),[],'','','','',Cmap, [], [], 1);
     set(gca,'FontSize',16); ylim([0,3]);
     ylabel('log_{10}(Fano Factor)','FontSize',16);
      xticklabel_rotate(linspace(.61,1.38,length(labs)),70,labs,'interpreter','none','FontSize',10);
      set(gcf,'color','w');
      
%   figure(2); clf;
%  barweb([ sna_COV; mouse_COV; yeast_COV; Mcell_COV ],...
%         [ sna_err; mouse_err; yeast_err; Mcell_err  ],[],...
%       '','','','CoV',Cmap);
%     ylim([.1,4]);
% 
%     figure(3); clf;
%  barweb([ sna_COV; mouse_COV; yeast_COV; Mcell_COV ],...
%         [ sna_err; mouse_err; yeast_err; Mcell_err  ],[],...
%       '','','','CoV',Cmap);
%     ylim([.08,.45]);
%     xlim([.6,.95]);

cmap2 = [cmap(1,:); cmap(2,:); cmap(4,:); cmap(6,:)];
figure(3); clf;  
clf; colordef black;
     barweb(all_COV(2:5),all_err(2:5),[],'','','','',cmap2);
     set(gca,'FontSize',15);
     ylabel('Cov','FontSize',15);   ylim([.08,.9]);
     set(gcf,'color','k');

   
 
 