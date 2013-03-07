function Image = LSM2image(finput,varargin)
%-------------------------------------------------------------------------
% LSM2image(Datas)
% LSM2image('LSMfile.lsm')
% LSM2image(Datas,z-section)
% LSM2image(Datas,z-section,position,fid)
% 
% reads the image from the variable Datas. 
%--------------------------------------------------------------------------

% default inputs
z = 1; 
p = 1;

% parse optional inputs
if nargin > 1
    z = varargin{1};
end
if nargin > 2
    p = varargin{2};
end
if nargin > 3
   fid = varargin{3};  
end
if nargin > 4
    disp('error, too many inputs');
end

% parse first input
if ischar(finput)
    Datas = LSMread(finput);
else
    Datas = finput;
end


%% Main code 

TIF = Datas.([ 'Stack' num2str(p)]).(['Image' num2str(z)]).TIF;
TIF.file = fid; 
Image = readLSMmatfile(TIF); 


