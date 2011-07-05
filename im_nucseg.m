
%%                  im_nucseg.m  Multi Channel
% Alistair Boettiger                                  Date Begund: 01/21/11
% Levine Lab, UC Berkeley                        Version Complete: 02/01/11 
% Functionally complete                             Last Modified: 02/07/11
% 
% 
%% Attribution:
% Feel free to use modify and distribute this code provided that you
% attribute Alistair Boettiger and Jacques Bothma for development and abide 
% by the provisions of the  Creative Commons License 3.0, BY-NC-SA.
% http://creativecommons.org/licenses/by-nc-sa/3.0/.
%
%
%
%%  Important Notes:
%  This version written for Mac.  To impliment in PC just change directory
% paths from using '/' to using '\'.  
% 
%  Before running, go scroll down to function setup and save the default
% parameters to your data folder
% 
%
%
%% Overview:
%  This code uses DNA staining to associate cytoplasmic domains with the
%  nearest nucleus.  High reslolution mRNA FISH localizes transcripts
%
%
%% Required subroutines
% fxn_nuc_seg.m  -- segmentation filter, identifies all nuclei
% fxn_nuc_reg.m -- expands nuclei to assign all regions of embryo to one
% nuclei or another.
% dotfinder.m -- locates dots using difference of gaussians and watershed
% DuplicateDots.m -- compares layers to ID duplicate dots
% vect2rast.m -- simple vector to raster conversion, called by DuplicateDots
% 
%% Updates: 
%  02/07/11 moved plotting 


function varargout = im_nucseg(varargin)
% IM_NUCSEG M-file for im_nucseg.fig
%      IM_NUCSEG, by itself, launches the GUI
%
%      IM_NUCSEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IM_NUCSEG.M with the given input arguments.
%
%      IM_NUCSEG('Property','Value',...) creates a new IM_NUCSEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before im_nucseg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to im_nucseg_OpeningFcn via varargin.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help im_nucseg
% Last Modified by GUIDE v2.5 10-Mar-2011 12:19:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @im_nucseg_OpeningFcn, ...
                   'gui_OutputFcn',  @im_nucseg_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before im_nucseg is made visible.
function im_nucseg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to im_nucseg (see VARARGIN)
   handles.output = hObject; % Choose default command line output for im_nucdots_v5
   
  % Some initial setup 
      % Folder to save .mat data files in for normal script function.  
     handles.fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';
     handles.step = 0;  % starting step is step 0 
     set(handles.stepnum,'String',handles.step); % change step label in GUI
     handles.output = hObject; % update handles object with new step number
     guidata(hObject, handles);  % update GUI data with new handles
     setup(hObject, eventdata, handles); % set up labels and default values for new step
     guidata(hObject, handles); % update GUI data with new labels        
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes im_nucseg wait for user response (see UIRESUME)
% uiwait(handles.figure1);




%=========================================================================%
%                          Primary Analysis Section                       %      
%=========================================================================%
% % All of the functional processing script is in this function
function run_Callback(hObject, eventdata, handles)
step = handles.step;

% Step 0: Load Data into script
if step == 0;
    disp('running...'); tic
    
     handles.NucChn = str2double(get(handles.in1,'String'));  % 
     handles.imsize = str2double(get(handles.in2,'String'));
     % save([handles.fdata,'/','test']);
     % load([handles.fdata,'/','test']);
 
    handles.output = hObject; % update handles object with new step number
   % guidata(hObject, handles);  % update GUI data with new handles
    handles = imload(hObject, eventdata,handles); % load new embryo
   guidata(hObject, handles);  % update GUI data with new handles
   
      
    
    
    toc
end

% Step 1: Nuclear Threshold
% uses fxn: fxn_nuc_seg
if step == 1;
% load appropriate data
    disp('running step 1...'); tic
    minN = str2double(get(handles.in1,'String'));  % 
    imblur = str2double(get(handles.in2,'String'));
    sigmaE = str2double(get(handles.in3,'String'));
    sigmaI = str2double(get(handles.in4,'String'));
    FiltSize = str2double(get(handles.in5,'String'));
  
    I = handles.Nucs; 
      
  % get threshold image 'bw' and nuclei centroids 'cent'  
    [handles.bw,handles.cent] = fxn_nuc_seg(I,minN,sigmaE,sigmaI,FiltSize,imblur);
   
 % Save data values  
 %      handles.output = hObject; guidata(hObject, handles);   
     guidata(hObject, handles);  % update GUI data with new handles 
    toc;
end
 
% Step 2: Get Region for each Nuclei
% uses fxn  fxn_nuc_reg
if step == 2;   
    tic
    disp('running step 2...');
    Mthink = str2double(get(handles.in1,'String'));  % 
    Mthin = str2double(get(handles.in2,'String'));
    Imnn = str2double(get(handles.in3,'String'));
    [NucLabeled,Nuc_overlay,conn_map,cell_bords] = fxn_nuc_reg(handles.Nucs,handles.bw,Mthink,Mthin,Imnn);  

    
    
    % save([handles.fdata,'/','test']);
    % load([handles.fdata,'/','test']);
    
    [h,w] = size(NucLabeled);
    
    Cell_bnd = false(h,w);
    Cell_bnd(cell_bords) = 1;
    Cell_bnd = bwareaopen( Cell_bnd,100);
    
    int = class(handles.Nucs);
    if strcmp(int,'uint16')
        Nucs = handles.Nucs + uint16(2^16*Cell_bnd);
    elseif strcmp(int,'uint8')
         Nucs = handles.Nucs + uint8(2^8*Cell_bnd);
    end
        
    
    figure(21); clf; imagesc(Nucs); colormap('hot');    
    
    
       
    handles.NucLabeled = NucLabeled; 
    handles.Cell_bnd = Cell_bnd; 
    handles.conn_map = conn_map; 
    guidata(hObject, handles);  % update GUI data with new handles  
    toc
end



 % export data
 if step == 3
 tic
    fout = get(handles.fout,'String');
    fname = get(handles.in1,'String');
    disp(['exporting data to ',fout,fname,'...']); 
    
        
    NucLabeled = handles.NucLabeled;
    nuc_cents = handles.cent; 
    Nucs = handles.Nucs; 
    conn_map = handles.conn_map; 
    Cell_bnd = handles.Cell_bnd;
    
    save([fout,'/',fname,'_nucdata.mat'],...
        'NucLabeled','nuc_cents','conn_map','Cell_bnd','Nucs'); 
  
    
    guidata(hObject, handles); 
 toc
   disp('data saved'); 
     
 end
 
 
%========================================================================%
 %  end of functional processing script
 % The rest of this code is GUI manipulations







% --- Executes on button press in VarButton.
function VarButton_Callback(hObject, eventdata, handles)



% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %
%                        File managining scripts                          %  
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %
% This function sets up the new steps with the appropriate input labels and
% defalut label parameters


  % handles.fdata = '/Users/alistair/Documents/Berkeley/Levine_Lab/ImageProcessing/';

function setup(hObject,eventdata,handles)
 if handles.step == 0; 
       load([handles.fdata, 'imnucseg_pars0']); 
       % pars = {'3','512',' ',' ',' ',' '}; save([handles.fdata,'imnucseg_pars0'], 'pars' );
        set(handles.in1label,'String','Nuclei channel');
        set(handles.in1,'String', pars(1));
        set(handles.in2label,'String','imsize');
        set(handles.in2,'String', pars(2));
       set(handles.in3label,'String',' ');
        set(handles.in3,'String', pars(3));
        set(handles.in4label,'String',' ');
        set(handles.in4,'String', pars(4));
        set(handles.in5label,'String',' ');
        set(handles.in5,'String', pars(5));
        set(handles.in6label,'String',' ');
        set(handles.in6,'String', pars(6));
            set(handles.VarButtonName,'String',''); 
        dir = {
       'Load max-projected image from folder.'} ;
        set(handles.directions,'String',dir); 
 end

 
 
   if handles.step == 1; 
       load([handles.fdata,'/','imnucseg_pars1']);
       %  pars = {'100','4','20','23','30',' '};  save([handles.fdata,'imnucseg_pars1'], 'pars' );
        set(handles.in1label,'String','min Nuc size'); % number of pixels in filter (linear dimension of a square)
        set(handles.in1,'String', pars{1});
        set(handles.in2label,'String','Imblur'); % width of Gaussian in pixels
        set(handles.in2,'String',pars{2});
        set(handles.in3label,'String','Excitation Width');
        set(handles.in3,'String',pars{3}); 
        set(handles.in4label,'String','Inhibition Width');
        set(handles.in4,'String', pars{4});
        set(handles.in5label,'String','FilterSize');
        set(handles.in5,'String', pars{5});
        set(handles.in6label,'String',' ');
        set(handles.in6,'String', pars{6});  
       dir = {
        'Step 1: Find nuclei.  Uses a difference of Gaussian filter with';
        'a min nucleus size filter and watershed splitter'}; 
        set(handles.directions,'String',dir);
  end      
  if handles.step == 2;  % nuclei segmentation
    load([handles.fdata,'/','imnucseg_pars2']); 
    % pars = {'45','3','2','','','',''};  save([handles.fdata,'imnucseg_pars2'], 'pars' );
        set(handles.in1label,'String','thicken nuclei'); 
        set(handles.in1,'String', pars{1});
        set(handles.in2label,'String','thin boundaries');
        set(handles.in2,'String', pars{2});
        set(handles.in3label,'String','erode'); 
        set(handles.in3,'String', pars{3});
        set(handles.in4label,'String',' ');
        set(handles.in4,'String', pars{4}); 
        set(handles.in5label,'String',' ');
        set(handles.in5,'String', pars{5}); 
        set(handles.in6label,'String',' ');
        set(handles.in6,'String', pars{6});  
                dir = {'Step 2: Map nuclear region';
    'nuclei expand until they collide.  Borders are assigned to different nuclei'} ;
        set(handles.directions,'String',dir); 
  end
  
  
  
  if handles.step == 3;
     load([handles.fdata,'/','imnucseg_pars3']); 
     % pars = {' ',' ',' ',' ',' ',' '}; save([handles.fdata,'imnucseg_pars3'], 'pars' );
     
     froot = get(handles.froot,'String');
     emb = get(handles.embin,'String');
     fname = [froot,'_',emb];  
     
     
        set(handles.in1label,'String','Save Name'); 
        set(handles.in1,'String', fname);
        set(handles.in2label,'String',' ');
        set(handles.in2,'String', pars{2});
        set(handles.in3label,'String',' '); 
        set(handles.in3,'String', pars{3}); 
        set(handles.in4label,'String',' ');
        set(handles.in4,'String', pars{4});
        set(handles.in5label,'String',' ');
        set(handles.in5,'String', pars{5});
        set(handles.in6label,'String',' ');
        set(handles.in6,'String', pars{6});
        set(handles.VarButtonName,'String','Manual Reg Select');
   dir = {'Step 3: Save Data'};
  set(handles.directions,'String',dir); 
  end
  
guidata(hObject, handles); % update GUI data with new labels





% --- Executes on button press in savePars.
function savePars_Callback(hObject, eventdata, handles)
   % record the values of the 6 input boxes for the step now showing
     p1 = get(handles.in1,'String');  
     p2 = get(handles.in2,'String');  
     p3 = get(handles.in3,'String');  
     p4 = get(handles.in4,'String');  
     p5 = get(handles.in5,'String');  
     p6 = get(handles.in6,'String');  
     
     try % for some reason the parameters are sometimes retrived as cells instead of strings
         % we need to make sure they are strings.  
       pars = {p1{:}, p2{:}, p3{:}, p4{:}, p5{:}, p6{:}}; % cell array of strings
     catch
        pars = {p1, p2, p3, p4, p5, p6}; % cell array of strings
     end 
  % Export parameters 
     stp_label = get(handles.stepnum,'String');     
     savelabel = ['imnucseg_pars',stp_label];  
     % labeled as nucdot_parsi.mat where "i" is the step number 
     save([handles.fdata, savelabel], 'pars');        % export values
     disp([handles.fdata, savelabel]);
     pars 

%      save([handles.fdata,'/','test']);
%      load([handles.fdata,'/','test']);
     
guidata(hObject, handles); % update GUI data with new labels

% ----------------------STEP CONTROLS----------------------- %
% interfaces to main analysis code 
% --- Executes on button press in nextstep.
function nextstep_Callback(hObject, eventdata, handles)
handles.step = handles.step + 1; % forward 1 step
 set(handles.stepnum,'String', handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
    setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels

% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
handles.step = handles.step-1; % go back a step
 set(handles.stepnum,'String',handles.step); % Change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles); % update GUI data with new handles
    setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels
% -------------------------------------------------------- %


% --- Executes on button press in LoadNext.
function LoadNext_Callback(hObject, eventdata, handles)
    embn = handles.emb + 1;  % update embryo number
    if embn<10
        emb = ['0',num2str(embn)];
    else
        emb = num2str(embn);
    end
    set(handles.embin,'String',emb); % update emb number field in GUI 
    handles.emb = emb; % update emb number in handles structure
   [handles] = imload(hObject, eventdata, handles); % load new embryo
   guidata(hObject, handles); % save for access by other functions

% Reset step to step 1; 
    handles.step = 1; % forward 1 step
    set(handles.stepnum,'String', handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
    setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels



        %========== change source images ================%
function [handles] = imload(hObject, eventdata,handles)

disp('loading image...');
 % handles.source
    handles.fin = get(handles.source,'String'); % folder
    handles.fname = get(handles.froot,'String'); % embryo name
    handles.embn = get(handles.embin,'String'); % embryo number
    handles.emb = str2double(handles.embn);
    
    % fname =  [handles.fin,'/',handles.fname,'_',handles.embn,'_max.tif']; 
    fname =  [handles.fin,'/','max_',handles.fname,'_',handles.embn,'.tif']; 
    handles.Im =  imread(fname);
    
    
    h = size(handles.Im,1);
    Nuc = handles.Im(:,:,handles.NucChn);
    scale = handles.imsize/h;
    handles.Nucs = imresize(Nuc,scale); 
    disp(['rescaling to ',num2str(scale*100,3), ' percent']); 
    
try
    Imini = imresize(handles.Im,scale);
    figure(1); clf; imshow(Imini);    
catch 
    try % try just nuclear layer.  
       figure(1); clf; imagesc(Nuc); colormap gray; 
    catch me
         disp(me.message);
    end
end
    
    
    
    handles.output = hObject; 
    guidata(hObject,handles);% pause(.1);
    disp('image loaded'); 
        %====================================================%



% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %



% Automatically return the program to step 0 if the image source directory,
% file name, or image number are changed.  

function froot_Callback(hObject, eventdata, handles)
 handles.step = 0;  % starting step is step 0 
     set(handles.stepnum,'String',handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
     setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels


function embin_Callback(hObject, eventdata, handles)
 handles.step = 0;  % starting step is step 0 
     set(handles.stepnum,'String',handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
     setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels



function source_Callback(hObject, eventdata, handles)
 handles.step = 0;  % starting step is step 0 
     set(handles.stepnum,'String',handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
     setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels


% Open file browser to select source folder 
function SourceBrowse_Callback(hObject, eventdata, handles)
 sourcefile = uigetdir; % prompts user to select directory
  set(handles.source,'String',sourcefile);





%% GUI Interface Setup
% The rest of this code just sets up the GUI interface


% --- Outputs from this function are returned to the command line.
function varargout = im_nucseg_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function embin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function source_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function froot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in2_Callback(hObject, eventdata, handles)
function in2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function in3_Callback(hObject, eventdata, handles)
function in3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in1_Callback(hObject, eventdata, handles)
function in1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function in5_Callback(hObject, eventdata, handles)
function in5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in6_Callback(hObject, eventdata, handles)
function in6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function in4_Callback(hObject, eventdata, handles)
function in4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function fout_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function fout_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

