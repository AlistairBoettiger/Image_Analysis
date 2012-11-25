function uiopen(type,direct)
% UIOPEN overloaded for custom Files. Do not change the file name of this
% file. Remember you are overloading uiopen inside toolbox/matlab/uitools
%

global myImage

%----Your file extension -----v
if ((~isempty(findstr(type,'.tif'))) && (direct))
    %-------------------------------------------------
    % Your function that will open/run this file type
    %-------------------------------------------------
    myImage = imread(type);
    %figure; imagesc(myImage);
    
    disp('Drag-Drop worked... ')
    %-------------------------------------------------
else
    %----------DO NOT CHANGE---------------------------
    presentPWD = pwd;
    cd([matlabroot '/toolbox/matlab/uitools']);
    strn = ['uiopen(''' type ''',' num2str(direct) ')'];
    eval(strn);
    cd(presentPWD);
    %----------DO NOT CHANGE---------------------------
end
%-------------------------------------------------------------------------