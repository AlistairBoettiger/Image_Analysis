% Imviewer
%
% Alistair Boettiger                                   Date Begun: 12/02/08
%                                                Version Complete: 01/19/09
% Levine Lab                                        Last Modified: 02/06/10
% 
%
% Takes raw data stacks from Leica microscope, projects z-stacks,
% superimposes colors, and the orients the image along the major axis.
% Enabled to autocycle through stack. 

% Updates:
% 02-06-10 updated to use however many sections are in stack if max slice
% is set to zero.  
%
%


function varargout = imviewer(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imviewer_OpeningFcn, ...
                   'gui_OutputFcn',  @imviewer_OutputFcn, ...
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


% --- Executes just before imviewer is made visible.
function imviewer_OpeningFcn(hObject, eventdata, handles, varargin)
handles.step = 1; 
setup(hObject,eventdata,handles); 
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = imviewer_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
step = handles.step; 
if step == 1;
    [hObject,handles] = step1(hObject, eventdata, handles);
end
if step == 2;
    [hObject,handles] = step2(hObject, eventdata, handles);
end
if step == 3;
     [hObject,handles] = step3(hObject, eventdata, handles);
end
 handles.output = hObject;
guidata(hObject, handles); 

%=========================================================================%
% STEP 1: PROJECT SELECTED RAW DATA FILES
function [hObject,handles] = step1(hObject, eventdata, handles) 
    % Input information
    fin = get(handles.fdir,'String');  % folder containing raw image data 
    handles.name =  get(handles.froot,'String'); %  file root 
    handles.imnum = get(handles.imnumber,'String'); %  image number as string

    % find labeled channels
    handles.chl = cell(3,1);
    handles.chl{3} = get(handles.bluename,'String');
    handles.chl{1} = get(handles.redname,'String');
    handles.chl{2} = get(handles.greenname,'String');
    
    chns = ['0'; '1'; '2'];
    if isempty(handles.chl{1});
        chns(1) = [];
    end
    if isempty(handles.chl{2});
        chns(2) = [];
    end
    if isempty(handles.chl{3});
        chns(3) = [];
    end

    save test
    % choose slices to use in quantification
    startslice = str2double(get(handles.in1,'String'));
    endslice = str2double(get(handles.in2,'String'));
    subset = [startslice,endslice]; 
    nuc_start = str2double(get(handles.in3,'String'));
    nuc_end = str2double(get(handles.in4,'String'));
    nuc_chn = str2double(get(handles.in5,'String'));
    nuc_dat = [nuc_start,nuc_end,nuc_chn]; 

    save test;
    Io = fxnproject_color(fin,handles.name,chns,subset,handles.imnum,nuc_dat); % project image from raw data
    handles.Io = Io;  save test;  % need to keep image
    handles.output = hObject;
    guidata(hObject, handles); 
    disp('image projected'); 
    
   

% STEP 2: CLEAN AND ORIENT IMAGE
    function [hObject,handles] = step2(hObject,eventdata,handles)
    % Image cleaning parameters
    emb_chn = str2double(get(handles.in1,'String'));
    objT = str2double(get(handles.in2,'String'));
    sigma= str2double(get(handles.in3,'String'));
    imdil = str2double(get(handles.in4,'String')); 
    minO= str2double(get(handles.in5,'String'));
    

    save test; 
    I = fxnclean_v2(handles.Io,emb_chn,objT,sigma,imdil,minO); % clean and orient image
   % figure(21); close; figure(20); close; % surpress fxn script figure output
    handles.I = I;
    figure(1); clf; imshow(I);
    guidata(hObject, handles); 
    handles.output = hObject;
    
 % STEP 3: SAVE DATA
   function [hObject,handles] = step3(hObject,eventdata,handles)
      [hObject,handles] = cupdate_Callback(hObject, eventdata, handles);
      outname = get(handles.outlabel,'String');
      fout = get(handles.foutput,'String'); 
      imwrite(handles.Ic,[fout,'/', outname,'_'...
               handles.imnum,'.tif'],'tif');
%=========================================================================%

function setup(hObject,eventdata,handles);  
if handles.step == 1; 
    load imviewer_1;
    % pars = [1,10,1,4,2]; save imviewer_1 pars;
       % set(handles.directions,'String','Step 1: Load and project raw image data');
        set(handles.in1label,'String','Start Slice');
        set(handles.in1,'String', num2str(pars(1)));
        set(handles.in2label,'String','End Slice');
        set(handles.in2,'String', num2str(pars(2)));
        set(handles.in3label,'String','Nuclei start slice');
        set(handles.in3,'String',num2str(pars(3)) );
        set(handles.in4label,'String','Nuclei end slice');
        set(handles.in4,'String',num2str(pars(4)) );
        set(handles.in5label,'String','Nuclei channel'); 
        set(handles.in5,'String',num2str(pars(5)) );
        set(handles.in5,'BackgroundColor','white'); 
        %set(handles.in3,'BackgroundColor',[.7 .7 .7]); 
        % set(handles.in4,'BackgroundColor',[.7 .7 .7]); 
end
if handles.step == 2;
   load imviewer_2;
   % pars = [3,0,30,50,300]; save imviewer_2 pars; 
        set(handles.in1label,'String','Object Channel'); 
        set(handles.in1,'String', num2str(pars(1)) );
        set(handles.in2label,'String','Object Threshold');
        set(handles.in2,'String', num2str(pars(2)) );
        set(handles.in3label,'String','Sigma for blur'); 
        set(handles.in3,'String', num2str(pars(3)) );
        set(handles.in4label,'String','Im Dilate');
        set(handles.in4,'String', num2str(pars(4)) ); 
        set(handles.in4,'BackgroundColor','white'); 
        set(handles.in5label,'String','Min Object Size'); 
        set(handles.in5,'String', num2str(pars(5)) );
       %set(handles.in5,'BackgroundColor',[.7 .7 .7]); 
end

% --- Executes on button press in Next.
function Next_Callback(hObject, eventdata, handles)
handles.step = handles.step + 1;
 set(handles.stepnum,'String',handles.step);
    handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles); 
 
% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
handles.step = handles.step - 1;
 set(handles.stepnum,'String',handles.step);
     handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles); 
 
    

% ---- Update colors
function [hObject,handles] = cupdate_Callback(hObject, eventdata, handles)
C(1,1) = str2double(get(handles.ch1R,'String'));
C(1,2) = str2double(get(handles.ch1G,'String'));
C(1,3) = str2double(get(handles.ch1B,'String'));
C(2,1) = str2double(get(handles.ch2R,'String'));
C(2,2) = str2double(get(handles.ch2G,'String'));
C(2,3) = str2double(get(handles.ch2B,'String'));
C(3,1) = str2double(get(handles.ch3R,'String'));
C(3,2) = str2double(get(handles.ch3G,'String'));
C(3,3) = str2double(get(handles.ch3B,'String'));

Ic = handles.I;
% use new structure to 
Ic(:,:,1) = C(1,1)*handles.I(:,:,1) + C(2,1)*handles.I(:,:,2) + C(3,1)*handles.I(:,:,3);
Ic(:,:,2) = C(1,2)*handles.I(:,:,1) + C(2,2)*handles.I(:,:,2) + C(3,2)*handles.I(:,:,3);
Ic(:,:,3) = C(1,3)*handles.I(:,:,1) + C(2,3)*handles.I(:,:,2) + C(3,3)*handles.I(:,:,3);


[w l] = size(Ic(:,:,1));
psat1 = sum(sum((Ic(:,:,1) == 255)))/(w*l);
psat2 = sum(sum((Ic(:,:,2) == 255)))/(w*l);
psat3 = sum(sum((Ic(:,:,3) == 255)))/(w*l);

meanI1 = mean(mean(Ic(:,:,1)));
meanI2 = mean(mean(Ic(:,:,2)));
meanI3 = mean(mean(Ic(:,:,3)));
[meanI1, meanI2, meanI3]

if psat1 > .005
    disp('Warning: red channel saturating');
end

if psat2 > .005
    disp('Warning: green channel saturating');
end

if psat3 > .005
    disp('Warning: blue channel  saturating');
end


            figure(1); clf; 
            imshow(Ic);
            handles.Ic = Ic; 
            handles.output = hObject; 
            guidata(hObject, handles);

% --- Executes on button press in flipH.
function flipH_Callback(hObject, eventdata, handles)
 I = handles.I;
 I(:,:,1) = fliplr(I(:,:,1));
 I(:,:,2) = fliplr(I(:,:,2));
 I(:,:,3) = fliplr(I(:,:,3));
 handles.I = I; 
 figure(1); clf; imshow(I); 
 guidata(hObject, handles);


% --- Executes on button press in flipV.
function flipV_Callback(hObject, eventdata, handles)
      I = handles.I;
 I(:,:,1) = flipud(I(:,:,1));
 I(:,:,2) = flipud(I(:,:,2));
 I(:,:,3) = flipud(I(:,:,3));
 handles.I = I; 
figure(1); clf; imshow(I); 
 guidata(hObject, handles);        
            
            
% --- Executes on button press in DelObj.
function DelObj_Callback(hObject, eventdata, handles)

        
xy = []; % initially the list of points is empty
n = 0; % initiate number of points

    dir1 = 'Left mouse button picks points.'; 
    dir2 = 'Right mouse button picks last point.';
   disp({dir1; dir2});
    % set(handles.directions,'String',{[dir1];[dir2]});

    figure(1); hold on; 
    
    but = 1;
    while but == 1
        [xi,yi,but] = ginput(1);
        plot(xi,yi,'wo')
        n = n+1;
        xy(:,n) = [xi;yi];
    end
    % Interpolate with a smooth (spline) curve and finer spacing.
    ts = linspace(1,n,20);
    xys = spline(1:n,xy,ts);

    % Plot the interpolated curve and export data
    bnd_x = [xys(1,:), xys(1,1)]; 
    bnd_y = [xys(2,:), xys(2,1)];
   
   figure(1); hold on;  
   plot(bnd_x,bnd_y,'m-');

   
   I = handles.I;
   [h,w] = size(I(:,:,1));
 DeleteMe = poly2mask(bnd_x,bnd_y,h,w);
 

   figure(2); clf; imshow(DeleteMe);
   
   I1 = I(:,:,1);  I1(DeleteMe) = 0; 
   I2 = I(:,:,2);  I2(DeleteMe) =0;
   I3 = I(:,:,3);  I3(DeleteMe) = 0;
   In(:,:,1) = I1;
   In(:,:,2) = I2;
   In(:,:,3) = I3;
   
   figure(2); clf; imshow(In);
   
 handles.I = In;
 handles.Io = In; 
figure(1); clf; imshow(I); 
 guidata(hObject, handles);      

            

% --- Autocycle
function cycle_Callback(hObject, eventdata, handles)
  h = 1; 
   while h==1;
       try
         % run script
          [hObject,handles] = step1(hObject, eventdata, handles);
          Next_Callback(hObject, eventdata, handles);
          [hObject,handles] = step2(hObject, eventdata, handles);
          Next_Callback(hObject, eventdata, handles);
          [hObject,handles] = step3(hObject, eventdata, handles);
                  
         % advance to next embryo
          imnum = str2double(get(handles.imnumber,'String'));
          imnum = imnum + 1; % advance to next image
          if imnum < 10; imnum = ['0',num2str(imnum)]; 
          else imnum = num2str(imnum); end;
          handles.imnum = imnum;
          set(handles.imnumber,'String',imnum); 
          
          % restart at step 1
            handles.step = 1;  
            set(handles.stepnum,'String',handles.step);
            setup(hObject, eventdata, handles); 
            handles.output = hObject; 
            guidata(hObject, handles);
          
          disp(['cycle: ',handles.imnum]);
        catch
          disp('autocycle complete');
          break
        end
  end


% --- update step default parameters 
function sdefaults_Callback(hObject, eventdata, handles)
stp_label = get(handles.stepnum,'String');
savelabel = ['imviewer_',stp_label];
pars = zeros(1,5); 
 pars(1) = str2double(get(handles.in1,'String'));
 pars(2) = str2double(get(handles.in2,'String')) ;
disp(savelabel); % save test;
try; pars(3) = str2double(get(handles.in3,'String')); catch; end;
try; pars(4) = str2double(get(handles.in4,'String')); catch; end;
try; pars(5) = str2double(get(handles.in5,'String')); catch; end;
 disp(num2str(pars));
save(savelabel,'pars'); 


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
function imnumber_Callback(hObject, eventdata, handles)
handles.step =  1;
 set(handles.stepnum,'String',handles.step);
    handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles);
  
function froot_Callback(hObject, eventdata, handles);
    froot = get(handles.froot,'String');
    set(handles.outlabel,'String',froot);
    guidata(hObject, handles);
        
        
        
function froot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end
    
function imnumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function fdir_Callback(hObject, eventdata, handles)
function fdir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end




function in1_Callback(hObject, eventdata, handles)
function in1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function in3_Callback(hObject, eventdata, handles)
function in3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function in4_Callback(hObject, eventdata, handles)
function in4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function in2_Callback(hObject, eventdata, handles)
function in2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function redname_Callback(hObject, eventdata, handles)
function redname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function greenname_Callback(hObject, eventdata, handles)
function greenname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function bluename_Callback(hObject, eventdata, handles)
function bluename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function ch1R_Callback(hObject, eventdata, handles)
function ch1R_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch1G_Callback(hObject, eventdata, handles)
function ch1G_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch1B_Callback(hObject, eventdata, handles)
function ch1B_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch3R_Callback(hObject, eventdata, handles)
function ch3R_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch3G_Callback(hObject, eventdata, handles)
function ch3G_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch3B_Callback(hObject, eventdata, handles)
function ch3B_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ch2R_Callback(hObject, eventdata, handles)
function ch2R_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch2G_Callback(hObject, eventdata, handles)
function ch2G_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch2B_Callback(hObject, eventdata, handles)
function ch2B_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in5_Callback(hObject, eventdata, handles)
% hObject    handle to in5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in5 as text
%        str2double(get(hObject,'String')) returns contents of in5 as a double


% --- Executes during object creation, after setting all properties.
function in5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function foutput_Callback(hObject, eventdata, handles)
% hObject    handle to foutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of foutput as text
%        str2double(get(hObject,'String')) returns contents of foutput as a double


% --- Executes during object creation, after setting all properties.
function foutput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outlabel_Callback(hObject, eventdata, handles)
% hObject    handle to outlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outlabel as text
%        str2double(get(hObject,'String')) returns contents of outlabel as a double


% --- Executes during object creation, after setting all properties.
function outlabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
