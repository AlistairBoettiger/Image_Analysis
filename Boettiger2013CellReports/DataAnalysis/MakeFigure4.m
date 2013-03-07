%%                  Comparison of snail dose effects
% Alistair Boettiger                                   Date Begun: 12/24/11
% Zhuang / Levine Lab                               Last Modified: 09/26/12

%% Produces Cell Reports Fig 4
%% publication figure as of 09/26/12
% This version paths setup to run on Monet.

%%
clear all;


 folder = 'C:\Users\Alistair\Documents\Projects\mRNA_counting\Data\';
addpath( 'C:\Users\Alistair\Documents\Projects\GenCode\'); 
% needs fxn_fit_sigmoid.m


fit_sigs = 1;%  optional fit sigmoids to curves
cluster_hist = 1;% optional plot cluster histograms
 F = 16; % fontsize for graphs

%'m199b', m109a
slides = {'MP05a','s142', 'wt','wt B','wt C','wt young 1','wt young 2',     'control BAC','control BAC B',    '4xsna','4xsna B',   '1xsna','1xsna B'};
groupnames = {'2x sna','2x BAC','4x sna','1x sna'};
wt =1; snaBAC = 2; FourX = 3; OneX = 4; 

nucclass =  'cc14'; %  'cc13meta';%  'cc13'; % choose nuc class to filter on
nucclass2 = 'unset';% 'cc13';%  'cc14l'; % optional additional class to add
nucclass3 = 'unset';% 'cc13meta';%  'cc14l'; % optional additional class to add
dat_label = cell(100,1); %


% Initialize variables and storage arrays
    G = length(groupnames);
    S = length(slides);
    E = 200; % max number of embryos per group
    olddataset = 0;  % total embryo counter resets eachtime dataset changes
    t = 0;% counter for total embryos
    
    
Nmolecules = 0; Ncells = 0;
Ncells_group = zeros(G,1);
Nembs_group = zeros(G,1);
   
  CMap = [0,0,1;
        0,.75,1;
       1,.2,.65;
      1,0,0;
       1,0,0];

%     figure(2); clf; % 
%     CMap = zeros(G,3); 
%     for g=1:G
%         CMap(g,:) = [g/G,0,1-g/G];
%     end

    cov = NaN*ones(E,G); % array to store CoV's for sna expressing ceills 
    aveon   = NaN*ones(E,G); % array to store mean snail concentrations
    if fit_sigs == 1
        hillcoefs = NaN*ones(E,G); 
        offsets = NaN*ones(E,G); 
        sigmax = NaN*ones(E,G); 
        bkd = NaN*ones(E,G); 
        fits = cell(E,G);  
        celldat = cell(E,G);  
    end
     
    xmin = 0; xmax = 1; 
    
    figure(2); clf; set(gcf,'color','w'); 
    
for s=1:S % loop through chosen slides
    switch slides{s}
        
        % 4x sna 
        case '4xsna'
            date = '2012-03-24_4xsna/';  fname ='4xsna';  ver = '';% '_v3';%
            sf = '';   chns = 1;   skip = 19;     dataset = FourX;   
            
        case '4xsna B'          
            date = '2012-03-24_4xsna/';  fname ='4xsna_d';  ver = '';% '_v3';%
            sf = '';   chns = 1;   skip = 19;     dataset = FourX;     
         
            % 1 x sna
        case '1xsna'
            date = '2011-12/'; fname = 'snaD'; ver = '';
            sf = ''; chns = 1; skip = 19; % use skip to remove non-1x guys
            dataset = OneX;
            % note, not all embryos in this collection are 1x 
            
        case '1xsna B'
         date = '2011-12/'; fname = 'snaD'; ver = '';
            sf = ''; chns = 1; skip = 19; % use skip to remove non-1x guys
            dataset = OneX;
            % note, not all embryos in this collection are 1x 
            
        %  WT data
    case 'wt'
        date = '2011-02-17/';  fname ='MP10_22C_sna_y_d';  ver = '_vN';% '_v3';%
        sf = 'MP10_22C/';   chns = 1;   skip = 9;     dataset = wt;     
        
    case 'wt B'
     date = '2011-02-17/';    fname = 'MP05_22C_sna_y_c';     ver = '_vN';% '_v2';
        sf = 'MP05_22C/'; chns = 1;        skip = [3,5,7];        dataset = wt;

    case 'wt C'
        date ='2011-06-20/';  fname ='s07_MP05Hz_22C_c';    ver = '_vN';% ''; % 
        sf = 's07_MP05Hz/';  chns = 1;        skip = 17;      dataset = wt;
      
    case 's142'
        date = '2011-12/';        fname = 's142_sna';        ver = '_v2'; 
        sf = ''; vout = '';        skip = [6,7,9,12];        dataset = wt;
      
    case 'MP05a'
        date = '2011-02-17/';         fname = 'MP05_22C_sna_y';        ver = '_vN'; 
        sf = 'MP05_22C/';  vout = '';        skip = 11;        dataset = wt;
    
    case 'wt young 1'
        date ='2011-05-22/';        fname = 's05_MP06Hz';         ver = '_vN';% '_v2'; %
        sf = 's05_MP06/';  vout = '';        chns = 1;        skip = 17;      dataset = wt;
      
    case 'wt young 2'
        date ='2011-05-22/';        fname ='s05_MP06Hz_b';          ver = '_v3';% '_v2'; %
        sf = 's05_MP06/';  dataset = wt;  chns = 1;    skip = 17;   vout = '';% '_o2'; 
                       
          
        % Control BAC
    case 'control BAC'
        date = '2011-04_and_earlier/';        fname ='MP12Hz_snaD_22C';         chns = 1;
        sf = 'MP12Hz/'; skip = 15;        ver = '';        vout = '';        dataset = snaBAC;
        
    case 'control BAC B'
        date = '2011-04_and_earlier/';         fname ='MP12Hz_snaD_22C_b'; 
        sf = 'MP12Hz/'; chns = 1;         skip = 1;      ver = '';    vout = '';  dataset = snaBAC; 
      
    case 'm199b'
            date = '2012-07-27_MLslides/';    fname =  'm199b_fewissues'; % 'm199b_issues'; %
              ver = ''; sf = ''; vout =''; chns =1; skip = [5,7]; % partial z-scans = 'issues'.
              dataset =snaBAC; 
            
    case 'm109a'
            date = '2012-08-26_mp07/';    fname = 'm109a';  ver = '_v2';
            sf = ''; vout =''; chns =1; skip = [10,11]; % partial z-scans;
                dataset =snaBAC; 
        
    
    end % end of listed slides
    
   
    load([folder,date,fname,ver,'_slidedata',vout,'.mat'],'data'); 
    Nembs =sum( logical( 1-cellfun(@isempty,data)) );
     try
        load([folder,date,fname,'_agedata.mat'],'nclass');
    catch er
        nclass =  cell(1, Nembs);
        disp(er.message); 
        for nn=1:Nembs
            load([folder,date,fname,'_',sprintf('%02d',nn),'_nucdata.mat']);
            nclass{nn} = AgeClass; 
            if strcmp(AgeClass,' ')
                disp({'Embryo Age data missing for:' ; [date,fname];
                    ['please enter it here and resave ', fname,'_',...
                    sprintf('%02d',nn),'_nucdata.mat']});
                error('missing age class');
                AgeClass = 'cc14';   save([folder,date,fname,'_',sprintf('%02d',nn),'_nucdata.mat'],'AgeClass','Cell_bnd','NucLabeled','Nucs','bw','conn_map','nuc_cents');
            end
        end
        disp(nclass);
    end
    
    
    
    

    if dataset ~= olddataset
        i = 0;
        olddataset = dataset; % reset olddata set
    end

    for e =1:18 % e = 9
        try        
           % i = i+1; % i actually keeps track of embryo number, s tracks slide
            
              % Filter embryo data based on nuc class 
              nc1 = logical(1-isempty(strfind(nclass{e},nucclass))); % contains string
               nc2 = strcmp(nucclass2, nclass{e}); % exact match
               nc3 = strcmp(nucclass3, nclass{e});  % exact match           
            if (nc1 || nc2 || nc3) && isempty(find(e==skip,1))  % only record data meeting condition
                disp('embryo found');
                i = i+1; % i actually keeps track of embryo number, s tracks slide
                t = t+1;
                if e<10; emb = ['0',num2str(e)]; else emb = num2str(e); end;
                dat_label{t} = [date,sf,'max_',fname,'_',emb,'.tif']; 
                
                % align gradients 
               if fit_sigs == 1 % Sigmoidal fitting routine
                  xdat = data{e}.Data_sort(:,1)/1000; % renaming data, convert to microns
                  sna = data{e}.Data_sort(:,2)-min(data{e}.mu(:,1)); % subtract background 
                  grad = sna; L = length(grad); st = floor(L/5);
                  % truncate data from right to get only sigmoids and not full steps  
                  grad = grad(st:end);  gx = xdat(st:end);
                 n = 3; % initial guesses for sigmoid fits
                 theta = 2*mean(gx);   % threshold
                 A = max(grad); % fxn max
                 b=min(grad);   % offset / background
                 [p,fit] = fxn_fit_sigmoid(gx',grad',[n,theta,A,b],'r');
                 hillcoefs(i,dataset) = p(1);
                 offsets(i,dataset) = abs(p(2));
                 sigmax(i,dataset) = p(3); 
                 bkd(i,dataset) = p(4); 
                 fits{i,dataset} = [gx'; fit]; 
                 celldat{i,dataset} = [xdat';sna'];
                 
                 xmin = min(xmin, -offsets(i,dataset)); % for global plotting
                 xmax = max(xmax, max(xdat)); % for global plotting
               end
                
                 
                % extract mean and std of ON cells 
                if fit_sigs == 2; % makes some bad calls
                    pk = sigmax(i,dataset); 
                else
                    pk = max(data{e}.mu(:,1) - min(data{e}.mu(:,1)));
                    sna = data{e}.Data_sort(:,2)-min(data{e}.mu(:,1)); % subtract background  
                end
                 
             
              
                
              % split on and off populations
                sna_on = sna(sna>1/2*pk); % separate out greater than half-max
                sna_off = sna(sna<1/2*pk);
              % record coefficient of variance among on cells
                aveon(i,dataset) = nanmean(sna_on); 
                cov(i,dataset) = nanstd(sna_on)/nanmean(sna_on);
                
              % plot cluster hists   
              if cluster_hist == 1
                figure(2);  hold on;
                L = length(sna_on);
                plot(.7*rand(1,L)+t*ones(1,L),sna_on,'.','color',CMap(dataset,:),'MarkerSize',5);
                L2 = length(sna_off);
                plot(.7*rand(1,L2)+t*ones(1,L2),sna_off,'.','color',CMap(dataset,:)*.2+[.8,.8,.8],'MarkerSize',5);
              
                Ncells = Ncells + length(sna_on) + length(sna_off);
                Nmolecules = Nmolecules + sum(sna); 
              end
              
              Ncells_group(dataset) = Ncells_group(dataset)  + length(sna_on) + length(sna_off);
              Nembs_group(dataset) = i; 
%               % average
%               if ave_curves == 1
%                   mRNA{t} = sna;
                  
            end
        catch er
            disp(er.message);
        end
    end

end
figure(2); ylim([0,300]); xlim([0,t+2]);
xlabel('embryo number','FontSize',F); ylabel('mRNA count','FontSize',F);
set(gca,'FontSize',F);

title([num2str(Nmolecules,'%10.0f') ' mRNAs from ',num2str(Ncells), ' cells']); 
%% Summary Plots
names = cell(length(groupnames),1);

% cdf plots of CoV expression
    figure(6); clf;
        cov(cov==0) = NaN;
        cmap = lines(length(groupnames));    
        for g=1:length(groupnames)
            N1 = sum(true-isnan(cov(:,g)));
            names{g} = [groupnames{g},' N=',num2str(N1)];
            cov1 = cov(logical(true-isnan(cov(:,g))),g);
            try
            [f,x]=ecdf(cov1); stairs(x,f,'color',CMap(g,:),'linewidth',3); hold on; 
            catch er; disp(er.message); end
        end
    legend(names,'Location','Best'); xlabel('CoV','FontSize',14);  
    set(gca,'FontSize',14); xlim([0,.45]);

% Boxplot of CoV
    figure(4); clf; boxplot(cov,'width',.9,'labels',names,'colors',CMap, 'plotstyle','compact','whisker',10); % ,'plotstyle','compact'
    set(gcf,'color','w'); set(gca,'FontSize',14); ylim([0,.42]);
    ylabel('CoV');
    p = ranksum(cov(logical(true-isnan(cov(:,1))),1),cov(logical(true - isnan(cov(:,2))),2));
    title(['CoV for ',nucclass, ' embryos.  p = ',num2str(p,2)]); 
    hand = findobj(gca,'linewidth',4); set(hand,'linewidth',20);
    whisk = findobj(gca,'-not','color','k','-and','linewidth',.5);
    set(whisk,'LineStyle','--','linewidth',2);
    

% boxplot of average expression
    aveon(aveon==0) = NaN;
    figure(5); clf; boxplot(aveon,'width',.9,'labels',names,'colors',CMap, 'plotstyle','compact','whisker',10); % ,'plotstyle','compact'
    set(gcf,'color','w'); set(gca,'FontSize',14); ylim([0,300]);
    ylabel('mRNA count');
    [p,h] = ranksum(aveon(logical(true-isnan(aveon(:,1))),1),aveon(logical(true - isnan(aveon(:,2))),2));
    title(['Ave Expression for ',nucclass, ' embryos.  p = ',num2str(p,2)]); 
    hand = findobj(gca,'linewidth',4); set(hand,'linewidth',20);
    whisk = findobj(gca,'-not','color','k','-and','linewidth',.5);
    set(whisk,'LineStyle','--','linewidth',2);
    


%% Spatial Plots

figure(21); clf; figure(22); clf;
pts = 1000; % number of points for interpolation
stp = 20; % interval at which to plot spread bars
xs = linspace(xmin,xmax-xmin,pts); % for interp. must be positive
xo = xs + xmin; % for plotting, 0 is located at snail 'boundary'

for g=1:G
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
        figure(22);  plot(xo,cave,'color',CMap(g,:),'linewidth',3); hold on;
        plot(xo,cave-.5*cstd,'--','color',CMap(g,:)); hold on;
        plot(xo,cave+.5*cstd,'--','color',CMap(g,:)); hold on;
              
        X = [xo,xo]; Y = [cave-.5*cstd,cave+.5*cstd];
        Y(isnan(Y)) = 0; % can't have NaN's in plot 
        
        verts = [X',Y']; % 
        faces = [linspace(1,pts,pts),fliplr(linspace(pts+1,2*pts,pts))];      
        
        figure(21); 
        p = patch('Faces',faces,'Vertices',verts,'FaceColor',CMap(g,:),'EdgeColor','none');
        hold on;  alpha .7;
         plot(xo,cave,'color',CMap(g,:),'linewidth',2); hold on;
        
end
figure(22);
   xlim([-4E1,4E1]); xlabel('distance from boundary, (um)','FontSize',F);
   ylim([0,1.2*max(aveon(:))]);  ylabel('mRNA count','FontSize',F); 
   set(gca,'FontSize',F); % legend(groupnames); 
   set(gcf,'color','w');
   
   figure(21); 
   xlim([-4E1,4E1]); xlabel('distance from boundary, (um)','FontSize',F);
   ylim([0,1.2*max(aveon(:))]);  ylabel('mRNA count','FontSize',F); 
   set(gca,'FontSize',F); % legend(groupnames); 
      set(gcf,'color','w');
   
   