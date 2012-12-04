
%%                  imviewer_lsm.m  Multi Channel
% Alistair Boettiger                                  Date Begund: 02/12/11
% Levine Lab, UC Berkeley                        Version Complete: 02/12/11 
% Functionally complete                             Last Modified: 06/01/11
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
%
%
%% Required subroutines
% jacquestiffread and dependencies
% 
%% Updates: 
% 03/10/11: changed export name to max_fname from fname_max.
% it's easier to work with filenames with the number at the end, not in the
% middle.  
% 06/01/11 overhauled to run in a single step, keep parameters in local
% directory, only import one layer at a time.  

function varargout = imviewer_lsm(varargin)
% IMVIEWER_LSM M-file for imviewer_lsm.fig
%      IMVIEWER_LSM, by itself, launches the GUI
%
%      IMVIEWER_LSM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMVIEWER_LSM.M with the given input
%      arguments.
%
%      IMVIEWER_LSM('Property','Value',...) creates a new IMVIEWER_LSM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imviewer_lsm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imviewer_lsm_OpeningFcn via varargin.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imviewer_lsm
% Last Modified by GUIDE v2.5 03-Dec-2012 12:43:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imviewer_lsm_OpeningFcn, ...
                   'gui_OutputFcn',  @imviewer_lsm_OutputFcn, ...
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


% --- Executes just before imviewer_lsm is made visible.
function imviewer_lsm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imviewer_lsm (see VARARGIN)
   handles.output = hObject; % Choose default command line output for im_nucdots_v5
   
  % Some initial setup 
      % Folder to save .mat data files in for normal script function.  
     handles.fdata = 'C:\Users\Alistair\Documents\Projects\Snail Patterning\Code_data/';
 
     
     handles.step = 1;  % starting step is step 0 
     set(handles.stepnum,'String',handles.step); % change step label in GUI
     handles.output = hObject; % update handles object with new step number
     guidata(hObject, handles);  % update GUI data with new handles
     setup(hObject, eventdata, handles); % set up labels and default values for new step
     
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imviewer_lsm wait for user response (see UIRESUME)
% uiwait(handles.figure1);




%=========================================================================%
%                          Primary Analysis Section                       %      
%=========================================================================%
% % All of the functional processing script is in this function
function run_Callback(hObject, eventdata, handles)
step = handles.step;

% Step 1: Max Project nuclear channel at 1024  1024 resoultion
if step == 1; 
    disp('running step 1...'); 
  
    handles.fname = get(handles.in1,'String'); 
    firstc =  get(handles.in2,'String'); 
    lastc = get(handles.in3,'String');
    handles.first = eval(firstc{:});
    handles.last = eval(lastc{:}); 
    
    handles.folderout = get(handles.fout,'String');   
    
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
    [handles] = projectNsave(hObject, eventdata, handles); % load new embryo
    guidata(hObject, handles);  
end



function [handles] = projectNsave(hObject, eventdata, handles)    
   
    % processes inputs, then call projection file
    froot = get(handles.froot,'String');  % .lsm file name
    oname = get(handles.oname,'String'); % output file name
    emb = get(handles.embin,'String'); % input file number
    emb_out = get(handles.in5,'String'); % output file number
    fout = get(handles.fout,'String');  % save folder
    fsource = get(handles.source,'String'); % source folder containing data
    nmax = str2double(get(handles.in1,'String')); % max noise parameter
    firstc =  get(handles.in2,'String'); % first frame for max project
    lastc = get(handles.in3,'String'); % last frame for max project  
    str_chns = get(handles.in4,'String'); % RGB/CMYK channel order
    
    handles.first = eval(firstc{:});
    handles.last = eval(lastc{:}); 
    
    [handles] = RunProjection(hObject,eventdata,handles,...
            froot,oname,emb,emb_out,fout,fsource,nmax,str_chns);
    guidata(hObject, handles); 

% this function is spliced out of ProjectNSave to enable automatic filename
% detection    
function [handles] = RunProjection(hObject,eventdata,handles,...
            froot,oname,emb,emb_out,fout,fsource,nmax,str_chns)    
    
    global TIF;
    
    if emb_out{:} == ' ';
        emb_out = emb;
    else
        emb_out = emb_out{:};
    end
    

    first = handles.first;
    last = handles.last;
    
    
    N = str2double(emb);  % embryo number

    if N == 1 % only need to do this once.  
         jacquestiffread([fsource,'/',froot,'.lsm']);
        %parselsm([fsource,'/',froot,'.lsm']);
    end
    
     % froot = 's02_MP01_Hz_22C'; fsource = '/Volumes/Data/Lab Data/Raw_Data/2011-05-22/';
   try
    load([fsource,'/',froot,'.mat'])  
   catch er
       disp(er.message);
        jacquestiffread([fsource,'/',froot,'.lsm']);
        % parselsm([fsource,'/',froot,'.lsm']); % new version
   end
    filetemp= fopen(Datas.filename,'r','l');
    
       
    % find if data is uint16 or something else; 
try
    inttype =  ['uint',num2str(Datas.Stack1.Image1.IMG.bits)];
        disp(['data is ',inttype]); 
     channels =  Datas.Stack1.Image1.TIF.SamplesPerPixel;  
     w = Datas.Stack1.Image1.IMG.width;
     h = Datas.Stack1.Image1.IMG.height;
        disp(['Data contains ', num2str(channels),' channels']);
catch er
    disp(er.message); 
    disp('data is not an image stack');
end
      
    if str_chns{:} == ' '
        chns = 1:channels;
    else
        chns = eval(str_chns{:});
        channels = length(chns); 
    end   
   
    Zs = Datas.LSM_info.DimensionZ;
    last(last == 0) = Zs; 
    Imax = zeros(h,w,channels,inttype); 
   tic
   disp('writing data...'); 
   
    for i=1:Zs  % i = 47
      %  disp(['embryo ', num2str(N), '   stack-position', num2str(i) ]);
       % shorthand    
    TIF = Datas.([ 'Stack' num2str(N)]).(['Image' num2str(i)]).TIF;
    IMG = Datas.([ 'Stack' num2str(N)]).(['Image' num2str(i)]).IMG;
    TIF.file=filetemp;
    Im_layer = zeros(h,w,channels,inttype); %   eval([inttype,'(zeros(h,w,channels));']);
   
        %read the image channels
        for cc = 1:channels % c=2
            offset = 0;    
            c = chns(cc);
            TIF.StripCnt = c;
            IMG.data{c} = read_planeT(offset, IMG.width, IMG.height, c,TIF); 
                 
          %   figure(1); clf; imagesc(IMG.data{c})
            
            
%             %  UINT12 FIX
%            if  strfind('uint16',inttype) % uint12 get called uint16 without having correct scaling
%             IMG.data{c} = makeuint(IMG.data{c},16); % for uint12
%            end
           
            % check for screwed up offset
            [h,w] = size(IMG.data{c});
            % look in middle of data set and see if it's noise or signal
            sdata =  IMG.data{c}( floor(h/2*.9):floor(h/2*1.1), floor(w/2*.9):floor(w/2*1.1)  ); 
            isnoise = std(double(sdata(:)));

% % trouble shooting                       
%  handles.fdata = 'C:\Users\Alistair\Documents\Projects\Snail Patterning\Code_data/';            
%      save([handles.fdata,'/','test']);
%    load([handles.fdata,'/','test']);
            
           if isnoise > nmax %  && c<3
               offset = 1;
               IMG.data{c} = read_planeT(offset, IMG.width, IMG.height,c,TIF); 
               sdata =  IMG.data{c}( floor(h/2*.9):floor(h/2*1.1), floor(w/2*.9):floor(w/2*1.1)  ); 
               fixed = std(double(sdata(:)));
               disp(['offset error found in chn ', num2str(c), ' layer ',num2str(i),...
                   '  std=',num2str(isnoise,5), '  now=',num2str(fixed,5)] );

               if fixed >  nmax  % 
                   if TIF.BitsPerSample(1) == 16
                      IMG.data{c} = uint16(zeros(h,w));
                   else
                      IMG.data{c} = uint8(zeros(h,w));
                   end
                     disp(['Fix failed for in chn ', num2str(c),...
                         ' layer ',num2str(i),'  skipping this image...'] );
               end 
           end          
           
             % Compute max project
            % not enough memory to do one shot max project, need to do this
            % progressively.  Fortunately max doesn't care (unlike ave). 
            if i>first(c)-1 && i<last(c)+1
                 Im_layer(:,:,cc) = IMG.data{c};    % Insert into multicolor single layer, only if it's in the selected range.  
                 
%                  figure(1); clf; imagesc(imresize(Im_layer(:,:,1),.2)); colorbar;
%                  title(['embryo ', num2str(N), '   stack-position', num2str(i) ]); 
%                  pause(.1);
                 
                Imax(:,:,cc) = max( cat(3,Imax(:,:,cc),Im_layer(:,:,cc)),[],3);       
            end
            if channels == 2
                Im_layer(:,:,3) = Im_layer(:,:,2);
                Im_layer(:,:,2) = zeros(h,w,1,inttype); % eval([inttype,'(zeros(h,w,1));']); 
            end
        end  % close loop over colors
        imwrite(Im_layer,[fout,'/',oname,'_',emb_out,'_z', num2str(i),'.tif'],'tif');
     
            if channels == 4
                ImRGB = zeros(h,w,3,inttype);
                ImRGB(:,:,1) = Im_layer(:,:,1)+Im_layer(:,:,4);
                ImRGB(:,:,2) = Im_layer(:,:,2);
                ImRGB(:,:,3) = Im_layer(:,:,3)+Im_layer(:,:,4);
                imwrite(ImRGB,[fout,'/',oname,'_',emb_out,'_RBG_z', num2str(i),'.tif'],'tif');
            end
        
    clear Im_layer %  TIF IMG;
    end
    

% Can't write a 2 channel tif, need to convert to a 3 channel version.  
            if channels == 2
                Imax(:,:,3) = Imax(:,:,2); 
                Imax(:,:,2) = eval([inttype,'(zeros(h,w,1))']); 
            end
    imwrite(Imax,[fout,'/','max_',oname,'_',emb_out,'.tif'],'tif');
        if channels == 4
                maxRGB = zeros(h,w,3,inttype);
                maxRGB(:,:,1) = Imax(:,:,1)+Imax(:,:,4);
                maxRGB(:,:,2) = Imax(:,:,2);
                maxRGB(:,:,3) = Imax(:,:,3)+Imax(:,:,4);
                imwrite(maxRGB,[fout,'/max_',oname,'_',emb_out,'_RBG', num2str(i),'.tif'],'tif');
        end
    
    
    disp(['data written for embryo', emb]); 
    guidata(hObject, handles);  % update GUI data with new handles
    clear Imax; 
    
   %  fclose(filetemp);
    toc


%========================================================================%
 %  end of functional processing script
 % The rest of this code is GUI manipulations


 
% --- Executes on button press in ExtractLSM.
function ExtractLSM_Callback(hObject, eventdata, handles)
% hObject    handle to ExtractLSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    fsource = get(handles.source,'String'); % source folder containing data
    lsmdata = dir([fsource,filesep,'*.lsm']);
    
    for n=1:length(lsmdata)
        froot = regexprep(lsmdata(n).name,'.lsm',''); % get root
        oname = froot; 

        emb = get(handles.embin,'String'); % input file number
        emb_out = get(handles.in5,'String'); % output file number
        fout = get(handles.fout,'String');  % save folder

        nmax = str2double(get(handles.in1,'String')); % max noise parameter
        firstc =  get(handles.in2,'String'); % first frame for max project
        lastc = get(handles.in3,'String'); % last frame for max project  
        str_chns = get(handles.in4,'String'); % RGB/CMYK channel order

        handles.first = eval(firstc{:});
        handles.last = eval(lastc{:}); 

        [handles] = RunProjection(hObject,eventdata,handles,...
                froot,oname,emb,emb_out,fout,fsource,nmax,str_chns);
        guidata(hObject, handles); 
    end
        
% --- Executes on button press in AutoCycle.
function AutoCycle_Callback(hObject, eventdata, handles)

    froot = get(handles.froot,'String'); % root name of image
    emb = str2double(get(handles.embin,'String'));  
    handles.folderout = get(handles.fout,'String');   % where to save images

     Ttot = tic ;
     
     err = 0; 
    try  
     while err == 0;
         handles.emb = emb; % needed for imload to get the right embyro 

         if emb < 10
             embin = ['0',num2str(emb)];
         else
            embin = num2str(emb);
         end
            handles.fname = [froot,'_',embin];  % needed for savename. 
            set(handles.embin,'String',embin); 

        disp(['running embryo ',embin,'...']);
 
        handles.output = hObject; % update handles object with new step number
        guidata(hObject, handles);  % update GUI data with new handles
        [handles] = projectNsave(hObject, eventdata, handles); % load new embryo
        guidata(hObject, handles); 

        emb = emb + 1;          
     end
         
    catch error
        disp(error.message); 
        disp('Export finished.' ); tout = toc(Ttot); 
        disp([num2str(tout/60,3),' minutes']); 
   end
     
    

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %
%                        File managining scripts                          %  
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %
% This function sets up the new steps with the appropriate input labels and
% defalut label parameters

function setup(hObject,eventdata,handles)
  if handles.step == 1; 
       %load([handles.fdata,'/','imviewer_lsm_pars1']); 
       load('imviewer_lsm_pars1');
       
     froot = get(handles.froot,'String');
     emb = get(handles.embin,'String');
            
       % pars = {'1.2E4','[1,1,1,1]','[0,0,0,0]',' ',' ',' '}; save([handles.fdata,'imviewer_lsm_pars1'], 'pars' );
       
       % pars = {'1.2E4','[1,1,1,1]','[0,0,0,0]',' ',' ',' '};   save(['imviewer_lsm_pars1'], 'pars' );
        set(handles.in1label,'String','Max noise');
        set(handles.in1,'String', pars(1));
        set(handles.in2label,'String','starting frames');
        set(handles.in2,'String', pars(2));
       set(handles.in3label,'String','end frames');
        set(handles.in3,'String', pars(3));
        set(handles.in4label,'String','channel order');
        set(handles.in4,'String', pars(4));
        set(handles.in5label,'String','Alt. emb #');
        set(handles.in5,'String', pars(5));
        set(handles.in6label,'String',' ');
        set(handles.in6,'String', pars(6));
        %    set(handles.VarButtonName,'String',''); 
        dir = {'Step 1: Export layer data as tifs and max project between chosen';
            'starting and ending frames.  Use 0 for last frame to use all data.';
            'Optional: change channel order e.g. [4,3,2,1]'} ;
        set(handles.directions,'String',dir); 
 end
 
 
 





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
     savelabel = ['imviewer_lsm_pars',stp_label];  
     % labeled as nucdot_parsi.mat where "i" is the step number 
    % save([handles.fdata, savelabel], 'pars');        % export values
    % disp([handles.fdata, savelabel]);
    save(savelabel, 'pars'); 
    disp(savelabel);
     pars 

%      save([handles.fdata,'/','test']);
%      load([handles.fdata,'/','test']);
     
guidata(hObject, handles); % update GUI data with new labels


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %



% Automatically return the program to step 0 if the image source directory,
% file name, or image number are changed.  

function froot_Callback(hObject, eventdata, handles)
 handles.step = 1;  % starting step is step 0 
     set(handles.stepnum,'String',handles.step); % change step label in GUI
     froot = get(handles.froot,'String');
     set(handles.oname,'String',froot);  % automatically set export name equal input name 
     
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
     setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels


function embin_Callback(hObject, eventdata, handles)
 handles.step = 1;  % starting step is step 0 
     set(handles.stepnum,'String',handles.step); % change step label in GUI
    % set(handles.embin,'String','01'); % return embyro counter to 1
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
     setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels



function source_Callback(hObject, eventdata, handles)
 handles.step = 1;  % starting step is step 0 
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
function varargout = imviewer_lsm_OutputFcn(hObject, eventdata, handles) 
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


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)


function fout_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function fout_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton14_Callback(hObject, eventdata, handles)


function oname_Callback(hObject, eventdata, handles)

function oname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


