function  Datas = LSMread(filename,varargin)



% global TIF TEMPDat
filenameout=([filename(1:end-3) 'mat']);

%% set defaults values :
TIF.SampleFormat     = 1;
TIF.SamplesPerPixel  = 1;
TIF.BOS              = 'ieee-le';          %byte order string

TEMPDat.counterr=0;
TEMPDat.pastposi=0;
TEMPDat.rat=0;
TEMPDat.std=0;

if  isempty(strfind(filename,'.'))
    filename = [filename,'.tif'];
end

TIF.file = fopen(filename,'r','l');
if TIF.file == -1
    stkname = strrep(filename, '.tif', '.stk');
    TIF.file = fopen(stkname,'r','l');
    if TIF.file == -1
        error(['File "',filename,'" not found.']);
    else
        filename = stkname;
    end
end


[~, m] = fileattrib(filename);
% obtain the full file path:
filename = m.Name;
% find the file size in bytes:
% m = dir(filename);
% filesize = m.bytes;


%% read header
% read byte order: II = little endian, MM = big endian
byte_order = fread(TIF.file, 2, '*char');
if ( strcmp(byte_order', 'II') )
    TIF.BOS = 'ieee-le';                                % Intel little-endian format
elseif ( strcmp(byte_order','MM') )
    TIF.BOS = 'ieee-be';
else
    error('This is not a TIFF file (no MM or II).');
end


%% ---- read in a number which identifies file as TIFF format
tiff_id = fread(TIF.file,1,'uint16', TIF.BOS);
if (tiff_id ~= 42)
    error('This is not a TIFF file (missing 42).');
end

%% ---- read the byte offset for the first image file directory (IFD)
TIF.img_pos = fread(TIF.file, 1, 'uint32', TIF.BOS);

%%%%%%%%%%%%%    Main processing code     %%%%%%%%%%%%%%%%
count=1;
while  TIF.img_pos ~= 0 
    TIF.BOS;
    clear IMG;
    IMG.filename = filename;
    % move in the file to the first IFD
    status = fseek(TIF.file, TIF.img_pos, -1);
    if status == -1
        error('invalid file offset (error on fseek)');
    end
    %disp(strcat('reading img at pos :',num2str(TIF.img_pos)));

    %read in the number of IFD entries
    num_entries = fread(TIF.file,1,'uint16=>uint16', TIF.BOS);
    %disp(strcat('num_entries =', num2str(num_entries)));
 
    
    % Read and process each IFD entry, essentially reading in the header,
    % image by image. LSM data is stored only in the header of the first
    % file.
    ENT=[];
    for i = 1:num_entries    
        file_pos  = ftell(TIF.file);  % save the current position in the file
        TIF.entry_tag = fread(TIF.file, 1, 'uint16=>uint16', TIF.BOS);   % read entry tag
        ENT = [TIF.entry_tag;ENT];  % read entry
        [entry, TIF, TEMPDat] = readIFDentry(TIF,TEMPDat);
        %disp(strcat('reading entry <',num2str(TIF.entry_tag),'>'));
        switch TIF.entry_tag
            case 254
                TIF.NewSubfiletype = entry.val;      
            case 256         % image width - number of column
                IMG.width          = entry.val;
                TIF.width=[IMG.width];
            case 257         % image height - number of row
                IMG.height         = entry.val;
                TIF.ImageLength    = entry.val;
            case 258         % BitsPerSample per sample
                TIF.BitsPerSample  = entry.val;
                TIF.BytesPerSample = TIF.BitsPerSample / 8;
                IMG.bits           = TIF.BitsPerSample(1);
                %fprintf('BitsPerSample %i %i %i\n', entry.val);
            case 259         % compression
                if ( entry.val ~= 1 )
                    error(['Compression format ', num2str(entry.val),' not supported.']);
                end
            case 262         % photometric interpretation
                TIF.PhotometricInterpretation = entry.val;
                if ( TIF.PhotometricInterpretation == 3 )
                    warning('tiffread2:LookUp', 'Ignoring TIFF look-up table');
                end
            case 273         % strip offset
                TIF.StripOffsets   = entry.val;
                TIF.StripNumber    = entry.cnt;
%                TIF.StripOffsetsList   = [TIF.StripOffsetsList;TIF.StripOffsets']
                %fprintf('StripNumber = %i, size(StripOffsets) = %i %i\n', TIF.StripNumber, size(TIF.StripOffsets));
            case 277         % sample_per pixel
                TIF.SamplesPerPixel  = entry.val;
                %fprintf('Color image: sample_per_pixel=%i\n',
                %TIF.SamplesPerPixel);
            case 279         % strip byte counts - number of bytes in each strip after any compressio
                TIF.StripByteCounts= entry.val;
            case 284         %planar configuration describe the order of RGB
                TIF.PlanarConfiguration = entry.val;
            case 34412       % Zeiss LSM data
                LSM_info           = entry.val;
                disp('Zeiss LSM data')
            otherwise
                fprintf( 'Ignored TIFF entry with tag %i (cnt %i)\n', TIF.entry_tag, entry.cnt);
        end
        % move to next IFD entry in the file
        status = fseek(TIF.file, file_pos+12, -1); 
        if status == -1
            error('invalid file offset (error on fseek)');
        end       
    end
    %%%%%% end of reading in individual image header %%%%%%%%

    
    
    %%%%%%%%% read the position of the next IFD address:       
    TIF.img_pos = fread(TIF.file, 1, 'uint32', TIF.BOS); 
    if IMG.width~=128
        Data.(['I' num2str(count)]).TIF=TIF;
        Data.(['I' num2str(count)]).IMG=IMG;
        count=count+1;
    end
end

Data.LSM_info=LSM_info;
numberofstacks = (count-1)/Data.LSM_info.DimensionZ;
for i =1:numberofstacks
    for j=1:Data.LSM_info.DimensionZ      
        Datas.(['Stack' num2str(i)]).(['Image' num2str(j)]).TIF=Data.(['I' num2str((i-1)*Data.LSM_info.DimensionZ +j)]).TIF;
        Datas.(['Stack' num2str(i)]).(['Image' num2str(j)]).IMG=Data.(['I' num2str((i-1)*Data.LSM_info.DimensionZ +j)]).IMG;
    end
end

Datas.LSM_info=LSM_info;
Datas.filename=filename;
Datas.TEMPDat = TIF.TEMPDat;

save(filenameout,'Datas');
disp('data saved'); 
fclose(TIF.file);
end
%%%%%%%%%%%% End of Main processing code   %%%%%%%%%%%%%%%%



%% ==================sub-functions that reads an IFD entry:===================


function [nbBytes, matlabType] = convertType(tiffType)
switch (tiffType)
    case 1
        nbBytes=1;
        matlabType='uint8';
    case 2
        nbBytes=1;
        matlabType='uchar';
    case 3
        nbBytes=2;
        matlabType='uint16';
    case 4
        nbBytes=4;
        matlabType='uint32';
    case 5
        nbBytes=8;
        matlabType='uint32';
    case 7
        nbBytes=1;
        matlabType='uchar';
    case 11
        nbBytes=4;
        matlabType='float32';
    case 12
        nbBytes=8;
        matlabType='float64';
    otherwise
        error('tiff type %i not supported', tiffType)
end
end

%% ==================sub-functions that reads an IFD entry:===================

function  [entry, TIF, TEMPDat] = readIFDentry(TIF,TEMPDat)

%global TIF TEMPDat;
entry.tiffType = fread(TIF.file, 1, 'uint16=>uint16', TIF.BOS);
entry.cnt      = fread(TIF.file, 1, 'uint32=>double', TIF.BOS); %%%edit
%disp(['tiffType =', num2str(entry.tiffType),', cnt = ',num2str(entry.cnt)]);
[ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);

if entry.nbBytes * entry.cnt > 4
    %next field contains an offset:
    offset = fread(TIF.file, 1, 'uint32', TIF.BOS);
    status = fseek(TIF.file, offset, -1);
    if status == -1
        error('invalid file offset (error on fseek)');
    end
end

if TIF.entry_tag == 33629   % metamorph 'rationals'
    entry.val = fread(TIF.file, 6*entry.cnt, entry.matlabType, TIF.BOS);
elseif TIF.entry_tag == 34412  %TIF_CZ_LSMINFO
    [entry.val,TIF] = readLSMinfo(TIF);
    disp('LSMinforead')        
else
    if entry.tiffType == 5
        entry.val = fread(TIF.file, 2*entry.cnt, entry.matlabType, TIF.BOS);
    elseif TIF.entry_tag == 273 
       entry.val = fread(TIF.file, entry.cnt, entry.matlabType , TIF.BOS);
       entryvaltemp=entry.val;
       for i=1:length(entry.val)
          TEMPDat.rat= [TEMPDat.rat; TEMPDat.pastposi/entry.val(i)];
          TEMPDat.std= [TEMPDat.std;std(TEMPDat.rat)];
           if TEMPDat.pastposi/entry.val(i)>20*std(TEMPDat.rat)
               TEMPDat.counterr=TEMPDat.counterr+1;    
           end
           entry.val(i)= entry.val(i)+TEMPDat.counterr*(2^32); 
           TEMPDat.pastposi=entryvaltemp(i) ;       
       end
       TEMPDat.entryval = entry.val;   % Added export for troubleshooting
       TIF.TEMPDat=TEMPDat;
    else
           entry.val = fread(TIF.file, entry.cnt, entry.matlabType, TIF.BOS);  
    end
end

if ( entry.tiffType == 2 );
    entry.val = char(entry.val');
end

end


%% ==============Full-parse of LSM info:

function [R,TIF] = readLSMinfo(TIF)

% Read the LSM info
% global TIF ;
R.MagicNumber              = sprintf('0x%09X',fread(TIF.file, 1, 'uint32', TIF.BOS));
R.StructureSize            = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionX               = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionY               = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionZ               = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionChannels        = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionTime            = fread(TIF.file, 1, 'int32', TIF.BOS);
R.IntensityDataType        = fread(TIF.file, 1, 'int32', TIF.BOS);
R.ThumbnailX               = fread(TIF.file, 1, 'int32', TIF.BOS);
R.ThumbnailY               = fread(TIF.file, 1, 'int32', TIF.BOS);
R.VoxelSizeX               = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeY               = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeZ               = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginX                  = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginY                  = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginZ                  = fread(TIF.file, 1, 'float64', TIF.BOS);
R.ScanType                 = fread(TIF.file, 1, 'uint16', TIF.BOS);
R.SpectralScan             = fread(TIF.file, 1, 'uint16', TIF.BOS);
R.DataType                 = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetVectorOverlay      = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetInputLut           = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetOutputLut          = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetChannelColors      = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.TimeInterval             = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OffsetChannelDataTypes   = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetScanInformation    = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetKsData             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetTimeStamps         = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetEventList          = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetRoi                = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetBleachRoi          = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetNextRecording      = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DisplayAspectX           = fread(TIF.file, 1, 'float64', TIF.BOS);
R.DisplayAspectY           = fread(TIF.file, 1, 'float64', TIF.BOS);
R.DisplayAspectZ           = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OffsetMeanOfRoisOverlay  = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetTopoIsolineOverlay = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetLinescanOverlay    = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.ToolbarFlags             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetChannelWavelength  = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetChannelFactors     = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.ObjectiveSphereCorrection= fread(TIF.file, 1, 'float64', TIF.BOS);
R.OffsetUnmixParameters    = fread(TIF.file, 1, 'uint32', TIF.BOS);

%Read Scan information 
if (R.OffsetScanInformation>0)
    status = fseek(TIF.file, R.OffsetScanInformation, -1);
    if status == -1
        error('error on fseek');
    end

    count=sqrt(-1);   
    num=1;  
    while count~=0
        count=real(count);
        A1=sprintf('0x%09X',fread(TIF.file, 1, 'uint32', TIF.BOS));
        A2=fread(TIF.file, 1, 'uint32', TIF.BOS);
        A3=fread(TIF.file, 1, 'uint32', TIF.BOS);
        if A2==2
        A4 = fread(TIF.file, A3, '*char', TIF.BOS)';
        elseif A2==4
        A4 = fread(TIF.file, round(A3/4), 'int32', TIF.BOS)';
        elseif A2==5   
        A4 = fread(TIF.file, round(A3/8), 'float64', TIF.BOS)';
        else
        A4 = fread(TIF.file, A3, '*char', TIF.BOS)';
        end
        R.AllData.(['EntryNumber' num2str(num)]).FieldTag=A1;
        R.AllData.(['EntryNumber' num2str(num)]).DataType=A2;
        R.AllData.(['EntryNumber' num2str(num)]).Data=A4;   
        if (A2+A3==0) && ~(strcmp(A1,'0x0FFFFFFFF'))
            count=count+1;
        elseif (A2+A3==0) && (strcmp(A1,'0x0FFFFFFFF'))
            count=count-1;            
        else
        end
       num=num+1; 
    end
end

%Summarising pertinent scan information
R.DataSummary=[];
R.DataSummary.Fluorophores={};
R.DataSummary.ExcitationWavelength={};
for i=1:num-1,
    if strcmp(R.AllData.(['EntryNumber' num2str(i)]).FieldTag,'0x010000004')     
        R.DataSummary.Objective=R.AllData.(['EntryNumber' num2str(i)]).Data;    
    elseif strcmp(R.AllData.(['EntryNumber' num2str(i)]).FieldTag,'0x010000036')
      try  
        time = datenum('30-Dec-1899', 'dd-mmm-yyyy') + R.AllData.(['EntryNumber' num2str(i)]).Data;
        R.DataSummary.Time = datestr(time, 'mmmm dd, yyyy HH:MM');
      catch
          R.DataSummary.Time = '';% 
      end
    elseif strcmp(R.AllData.(['EntryNumber' num2str(i)]).FieldTag,'0x070000026')       
        R.DataSummary.Fluorophores(end+1,1)={R.AllData.(['EntryNumber' num2str(i)]).Data};
        disp(R.AllData.(['EntryNumber' num2str(i)]).Data)      
    elseif strcmp(R.AllData.(['EntryNumber' num2str(i)]).FieldTag,'0x090000001')
        
        R.DataSummary.ExcitationWavelength(end+1,1)={R.AllData.(['EntryNumber' num2str(i)]).Data};      
    else
    end  
end


end


