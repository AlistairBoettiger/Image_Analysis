%%                          imclean.m
% Alistair Boettiger                           Date Begun: 01/25/2012
% Zhuang Lab

% Stripped down version of previous imviewer_fast.fig

function varargout = imclean(varargin)
% IMCLEAN MATLAB code for imclean.fig
%      IMCLEAN, by itself, creates a new IMCLEAN or raises the existing
%      singleton*.
%
%      H = IMCLEAN returns the handle to a new IMCLEAN or the handle to
%      the existing singleton*.
%
%      IMCLEAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMCLEAN.M with the given input arguments.
%
%      IMCLEAN('Property','Value',...) creates a new IMCLEAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imclean_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imclean_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imclean

% Last Modified by GUIDE v2.5 01-Feb-2012 15:53:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imclean_OpeningFcn, ...
                   'gui_OutputFcn',  @imclean_OutputFcn, ...
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


% --- Executes just before imclean is made visible.
function imclean_OpeningFcn(hObject, eventdata, handles, varargin)
handles.step = 1; 
setup(hObject,eventdata,handles); 
set(handles.stepnum,'String','1');
handles.output = hObject;
guidata(hObject, handles);
 
A = zeros(100,200);
axes(handles.axes1); imagesc(A); axis off; colormap gray;
axes(handles.axes2); imagesc(A); axis off; colormap gray;
axes(handles.axes3); imagesc(A); axis off; colormap gray;
axes(handles.axes4); imagesc(A); axis off; colormap gray;

C = [1,0,0;
    0,1,0;
    0,0,1;
    0,0,0];
    set(handles.ctable,'Data',C);
    set(handles.ctable, 'ColumnEditable', [true true true]);
    handles.C = C;

% Choose default command line output for imclean
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imclean wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imclean_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes on button press in auto.
function auto_Callback(hObject, eventdata, handles)
  h = 1; 
   while h==1;
       try
         % run script
          [hObject,handles] = step1(hObject, eventdata, handles); % load image 
          NextStep_Callback(hObject, eventdata, handles);
          [hObject,handles] = step2(hObject, eventdata, handles); % clean and orient image
          NextStep_Callback(hObject, eventdata, handles);
          [hObject,handles] = step3(hObject, eventdata, handles); % save image 
                   
         % advance to next embryo
          imn = str2double(get(handles.imnum,'String'));
          imn = imn + 1; % advance to next image
          if imn < 10; imn = ['0',num2str(imn)]; 
          else imn = num2str(imn); end;
          handles.imn = imn;
          set(handles.imnum,'String',imn); 
          
          % restart at step 1
            handles.step = 1;  
            set(handles.stepnum,'String',handles.step);
            setup(hObject, eventdata, handles); 
            handles.output = hObject; 
            guidata(hObject, handles);
          
          disp(['cycle: ',handles.imn]);
       catch er
          disp('autocycle complete');
          break
        end
  end



% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)

step = handles.step; 
if step == 1;  % load image 
    [hObject,handles] = step1(hObject, eventdata, handles); 
    figure(20); clf; imagesc(handles.Io);
end
if step == 2; % clean and orient image
    [hObject,handles] = step2(hObject, eventdata, handles);
    figure(20); clf; imagesc(handles.I);
end
if step == 3; % save image
     [hObject,handles] = step3(hObject, eventdata, handles);
end
 handles.output = hObject;
guidata(hObject, handles); 



%%
% STEP 1: LOAD RAW IMAGE
    function [hObject,handles] = step1(hObject,eventdata,handles)
        fsource = get(handles.sourcef,'String');
        froot = get(handles.froot,'String');
        emb = get(handles.imnum,'String'); 
        filein = strcat(fsource,'/',froot,'_',emb,'.tif');
     %   save test;
        handles.Io = imread(filein);
      %  figure(20); clf; imagesc(handles.Io); % original image  
            
        guidata(hObject, handles); 
        handles.output = hObject;
        
        
% STEP 2: CLEAN AND ORIENT IMAGE
    function [hObject,handles] = step2(hObject,eventdata,handles)
    % Image cleaning parameters
    emb_chn = str2double(get(handles.in1,'String'));
    objT = str2double(get(handles.in2,'String'));
    objC = str2double(get(handles.in3,'String'));
    minO = str2double(get(handles.in4,'String'));

   % save test; 
  %~~~~~ 
       I = handles.Io;  
       inttype = class(I); inttype = str2double(inttype(5:end)) ;
      % locate object
      % objT = .15; objC = 12;
        bw = im2bw(I(:,:,emb_chn),objT); % Threshold image
        bw = imfill(bw,'holes'); % make objects solid        
        bw = bwareaopen(bw,minO); % remove objects of less than 100 pixels  
       % figure(2); clf; imagesc(1-bw); colorbar; colormap hot;
 
      sc = 2^4;  % should be parameter
      
      bw2 = imresize(bw,1/sc); 
      D = -bwdist(~bw2); % figure(2); clf; imagesc(D); colorbar; colormap hot;
      D = imfilter(D,fspecial('gaussian',30,5),'replicate');
      L = watershed(D); %  figure(2); clf; imagesc(L); colorbar; colormap hot;
      BW = bw2; BW(L==0)=0; 
      bw = imresize(BW,sc);
      % figure(2); clf; imshow(bw);

      % choose largest object, remove minors, then orient
         L = bwlabel(bw,8);   % create label matrix of objects in image
         
         imdata = regionprops(L,'Area', 'Orientation'); % measure properties of objects
         keep = find([imdata.Area]==max([imdata.Area])); % keep largest object
         theta = [imdata.Orientation]; % orientation of objects
         [w, l] = size(L);
         emb = zeros(w,l,3);
         obj_loc = ismember(L,keep); 
         obj_loc = imdilate(obj_loc,strel('disk',objC));
        %  figure(1); clf; imagesc(obj_loc);
         for k=1:3
            emb(:,:,k) = obj_loc; 
         end
         
      %   save test;
         I2 = makeuint(double(I).*double(emb),inttype); % deletes all other objects by mult zeros
       %  figure(21); clf; imagesc(I2); 
 
         I2 = imrotate(I2,360-theta(keep),'bilinear'); % rotate object        
     
      % resize image to zoom in on the object
         obj_loc = imrotate(obj_loc,360-theta(keep),'bilinear'); % align
         obj_loc = imdilate(obj_loc,strel('disk',10));  % add a little buffer
         L = bwlabel(obj_loc,8); 
         imdata = regionprops(L,'BoundingBox'); % compute boundaries
         c = round([imdata.BoundingBox]);
         try
         I3 = I2(c(2):c(4)+c(2),c(1):c(3)+c(1),:);    
         catch
             I3 = I2(c(2):c(4)+c(2)-1,c(1):c(3)+c(1)-1,:);
         end
         I = I3; % change output name
       %   figure(20); clf; imshow(I); title('cleaned,oriented embryo');
       
       try
axes(handles.axes1); imagesc(I(:,:,1)); axis off; colormap gray;
axes(handles.axes2); imagesc(I(:,:,2)); axis off; colormap gray;
axes(handles.axes3); imagesc(I(:,:,3)); axis off; colormap gray;
axes(handles.axes4); imagesc(I(:,:,4)); axis off; colormap gray;
       catch
       end
     
 %~~~~~~   
    handles.I_chn_ordered = I;  % memorize channel order to change color 
    handles.I = I; % updated field for output
    guidata(hObject, handles); 
    handles.output = hObject;
    
 % STEP 3: SAVE DATA
   function [hObject,handles] = step3(hObject,eventdata,handles)
      outname = get(handles.oroot,'String');
      fout = get(handles.outputf,'String'); 
      shiftnumber = 0; % str2double(get(handles.in1,'String')); 
      emb_number = str2double(get(handles.imnum,'String')); 
      
      emb = shiftnumber + emb_number;
      if emb < 10
          emb = ['0',num2str(emb)];
      else
          emb = num2str(emb);
      end
      disp('saving image...');
      
     %  save test;
      nameout = strcat(fout,'/', outname,'_',emb,'.tif');
      imwrite(handles.I,nameout,'tif');
      disp('image saved');
          guidata(hObject, handles); 
         handles.output = hObject;   
         
           

 %=========================================================================%

 
% --- Executes on button press in FlipH.
function FlipH_Callback(hObject, eventdata, handles)
 I = handles.I;
 I(:,:,1) = fliplr(I(:,:,1));
 I(:,:,2) = fliplr(I(:,:,2));
 I(:,:,3) = fliplr(I(:,:,3));
 handles.I = I; 
 figure(1); clf; imshow(I); 
 guidata(hObject, handles);


% --- Executes on button press in FlipV.
function FlipV_Callback(hObject, eventdata, handles)
      I = handles.I;
 I(:,:,1) = flipud(I(:,:,1));
 I(:,:,2) = flipud(I(:,:,2));
 I(:,:,3) = flipud(I(:,:,3));
 handles.I = I; 
figure(1); clf; imshow(I); 
 guidata(hObject, handles);    


% ---- Update colors
% --- Executes on button press in Recolor.
function Recolor_Callback(hObject, eventdata, handles)
    

    
 guidata(hObject, handles); 
 [hObject,handles] = cupdate_Callback(hObject, eventdata, handles);
  guidata(hObject, handles); 
 
 

% --- Executes when entered data in editable cell(s) in ctable.
function ctable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ctable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


function [hObject,handles] = cupdate_Callback(hObject, eventdata, handles)

 C = get(handles.ctable,'Data'); 
    
Ic = handles.I_chn_ordered;
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
            handles.I = Ic; % update output image with new colors
            handles.output = hObject; 
            guidata(hObject, handles);

 





function setup(hObject,eventdata,handles)
if handles.step == 1; 
    load imclean_1;
    % pars = [0,0,0,0]; save imclean_1 pars;
       % set(handles.directions,'String','Step 1: Load and project raw image data');
        set(handles.in1label,'String','');
        set(handles.in1,'String', num2str(pars(1)));
        set(handles.in2label,'String','');
        set(handles.in2,'String', num2str(pars(2)));
        set(handles.in3label,'String','');
        set(handles.in3,'String',num2str(pars(3)) );
        set(handles.in4label,'String','');
        set(handles.in4,'String',num2str(pars(4)) );

end
if handles.step == 2;
   load imclean_2;
   % pars = [1,.1,30,100]; save imclean_2 pars; 
        set(handles.in1label,'String','Object Channel'); 
        set(handles.in1,'String', num2str(pars(1)) );
        set(handles.in2label,'String','Object Threshold');
        set(handles.in2,'String', num2str(pars(2)) );
        set(handles.in3label,'String','Strel Close Object'); 
        set(handles.in3,'String', num2str(pars(3)) );
        set(handles.in4label,'String','Min Object Size');
        set(handles.in4,'String', num2str(pars(4)) ); 
        set(handles.in4,'BackgroundColor','white'); 
end

if handles.step == 3;
   load imclean_3;
   % pars = [0,0,0,0]; save imclean_3 pars; 
        set(handles.in1label,'String','Renumber from'); 
        set(handles.in1,'String', num2str(pars(1)) );
        set(handles.in2label,'String',' ');
        set(handles.in2,'String', num2str(pars(2)) );
        set(handles.in3label,'String',' '); 
        set(handles.in3,'String', num2str(pars(3)) );
        set(handles.in4label,'String',' ');
        set(handles.in4,'String', num2str(pars(4)) ); 
end

    guidata(hObject, handles); 
    handles.output = hObject;
%%
% --- Executes on button press in savepars.
function savepars_Callback(hObject, eventdata, handles)
stp_label = get(handles.stepnum,'String');
savelabel = ['imclean_',stp_label];
pars = zeros(1,4); 
 pars(1) = str2double(get(handles.in1,'String'));
 pars(2) = str2double(get(handles.in2,'String')) ;
pars(3) = str2double(get(handles.in3,'String'));
pars(4) = str2double(get(handles.in4,'String'));
 disp(num2str(pars));
save(savelabel,'pars'); 

%%

% --- Executes on button press in NextStep.
function NextStep_Callback(hObject, eventdata, handles)
handles.step = handles.step + 1;
 set(handles.stepnum,'String',handles.step);
    handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles); 


% --- Executes on button press in BackStep.
function BackStep_Callback(hObject, eventdata, handles)
handles.step = handles.step - 1;
 set(handles.stepnum,'String',handles.step);
     handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles); 
 


%% Automatically return to image 1, step 1 if these fields are changed

function imnum_Callback(hObject, eventdata, handles)
    handles.step =  1;
    set(handles.stepnum,'String',handles.step);
    handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles);
    
function sourcef_Callback(hObject, eventdata, handles)
   handles.step =  1;
   set(handles.stepnum,'String',handles.step);
   set(handles.imnum,'String','01');
   fin = get(handles.sourcef,'String');
   set(handles.outputf,'String',fin);
    handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles);
 
    function froot_Callback(hObject, eventdata, handles)
   handles.step =  1;
   set(handles.stepnum,'String',handles.step);
   set(handles.imnum,'String','01');
   rootin = get(handles.froot,'String');
   set(handles.oroot,'String',rootin);
    handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles);
  
%% More GUI Setup scripts, no need to edit
function outputf_Callback(hObject, eventdata, handles)

function oroot_Callback(hObject, eventdata, handles)

function oroot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function froot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
function imnum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    
    
function sourcef_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function outputf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function in1_Callback(hObject, eventdata, handles)

function in1_CreateFcn(hObject, eventdata, handles)

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

function in4_Callback(hObject, eventdata, handles)

function in4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function ctable_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ctable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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
