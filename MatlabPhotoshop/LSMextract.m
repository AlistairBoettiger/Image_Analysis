function LSMextract(LSMfile,varargin)    
%--------------------------------------------------------------------------
%  LSMextract(LSMfile)
% Extract tifs, tif Z-stacks, and/or z-stacks with multiple positions from
% an LSM file.  Saves the tifs in the same folder or a folder of your
% choice. 
% Compatiable with .lsm formats from Zeiss Zen Black 2012.  Earlier
% versions of Zen .lsm files also work but may produce errors for large
% files sizes (>10GB) or large numbers of positions.  
% 
%--------------------------------------------------------------------------
% Required Input:  LSMfile / string
%           - Full path and filename for the LSMfile to extract tiffs from.
% 
%--------------------------------------------------------------------------
% Optional Input: ('OptionName' / type / default)
% savename / string / same as original
%           - save name root for Tif files
% pathout / string / same as path in
%           - Directory to save tifs
% PositionList / double (vector)
%           - integer list of position in LSM file to extract
% PositionOffset / double 
%           - allows user to start numbering positions at a chosen value
%           (for example to have the output pickup from another experiment)
% maxZend / double / last Z-position
%           - vector of length channels.  Where each channel
%           should start or end.  
% maxZstart / double / first Z-position
%           - vector of length channels.  Where each channel
%           should start or end. 
% chns / double / vector
%           - which channels should be displayed and how they should be
%           ordered.  [3,2,1] flips the first and last channel (so that
%           what was blue in RGB is now red and v.v.). 

% Defaults for Optional Inputs
[pathin,savename] = extractpath(LSMfile);
pathout = pathin;
savename = regexprep(savename,'.lsm','');
maxZend = [];
maxZstart = [];
chns = [];
position_offset = 0; 
PositionList = [];

%--------------------------------------------------------------------------
%% Parse Variable Input Parameters
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'savename'
                savename = CheckParameter(parameterValue, 'string', 'savename');
            case 'pathout'
                pathout = CheckParameter(parameterValue, 'string', 'pathout');
            case 'PositionList'
                PositionList = CheckParameter(parameterValue,'nonnegative','PositionList'); 
            case 'PositionOffset'
                position_offset= CheckParameter(parameterValue, 'nonnegative', 'PositionOffset');
            case 'maxZend'
                maxZend= CheckParameter(parameterValue, 'nonnegative', 'maxZend');
            case 'maxZstart'
                maxZstart= CheckParameter(parameterValue, 'nonnegative', 'maxZstart');
            case 'chns'
                chns= parameterValue;
            otherwise
                error(['The parameter ''', parameterName,...
                    ''' is not recognized by the function, ''',...
                    mfilename '''.' '  See help ' mfilename]);
        end
    end
end


%%

% LSMfile = 'D:\Data\2013-01-14_no-prox\MP05_a.lsm'
Datas  = LSMread(LSMfile);
fid = fopen(Datas.filename,'r','l');

% Extract number of channels, number of positions, and number of Zs.
    NumPositions = length(fieldnames(Datas))-3;   
    inttype =  ['uint',num2str(Datas.Stack1.Image1.IMG.bits)];
        disp(['data is ',inttype]); 
     channels =  Datas.Stack1.Image1.TIF.SamplesPerPixel;  
        disp(['Data contains ', num2str(channels),' channels']);
    Zs = Datas.LSM_info.DimensionZ;

% Parse start and end for each z-stack appropriately
    if isempty(maxZend)
        maxZend = Zs*ones(1,channels);
    end
    if isempty(maxZstart)
        maxZstart = ones(1,channels);
    end

% if user does not specify channels, use all channels. 
    if isempty(chns) 
        chns = 1:channels;
    end   
    if channels < 4
        chnout = 3;
    else
        chnout = 4;
    end
   
    if isempty(PositionList)
        PositionList = 1:NumPositions;
    end
        

    
%% 
for p=1:PositionList
    disp(['reading data for position ',num2str(p),'...']);     
    tic;
    Imax = [];
    for z=1:Zs  % i = 24   
        Imtemp = LSM2image(Datas,z,p,fid); 
        if isempty(Imax)
            [h,w,Cs] = size(Imtemp);
            Imax = zeros(h,w,Cs,inttype);
        end
        for c = 1:channels % c=2  c=3
            if z>maxZstart(c)-1 && z<maxZend(c)+1                 
                Imax(:,:,c) = max(cat(3,Imax(:,:,c),Imtemp(:,:,c)),[],3);       
            end
        end  % close loop over colors
        
        % Save images 
        
        posNumberOut = sprintf('%02d',p+position_offset);
        Im_layer = zeros(h,w,chnout,inttype);    
        Im_layer(:,:,chns) = Imtemp; 
            imwrite(Im_layer,[pathout,'/',savename,'_',posNumberOut,'_z', sprintf('%03d',z),'.tif']);
            if chnout == 4
                ImRGB = Ncolor(Im_layer);
                imwrite(ImRGB,[pathout,'/',savename,'_',posNumberOut,'_RBG_z',sprintf('%03d',z),'.tif']);
            end     
    end % close loop over z-stack indices

    if Zs > 1 
     ImaxOut = zeros(h,w,chnout,inttype);    
     ImaxOut(:,:,chns) = Imax; 
     imwrite(ImaxOut,[pathout,'/','max_',savename,'_',posNumberOut,'.tif']);
      if channels == 4
        maxRGB = Ncolor(ImaxOut);
        imwrite(maxRGB,[pathout,'/max_',savename,'_',posNumberOut,'_RBG','.tif']);
      end
    end
  
    toc
    disp(['data written for position', num2str(p)]); 
end % close loop over positions
 fclose(fid);

