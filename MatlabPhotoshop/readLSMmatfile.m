function image = readLSMmatfile(TIF,varargin)
% 
%--------------------------------------------------------------------------
% Inputs
% TIF -- matlab structure defined in LSMread
% 
% Optional Inputs
% channels - vector indicating which channel(s) to extract.  
% 

% plane_nb % z-plane to read
Nchannels = TIF.SamplesPerPixel;
width = TIF.width;
height = TIF.ImageLength;
if nargin > 1
    channels = varargin{1};
else
    channels =1:Nchannels;
end


%determine the type needed to store the pixel values:
% all channels necessarily have the same bit-depth.  Don't need separate
% channel indices.  
switch(TIF.SampleFormat)
    case 1
        classname = sprintf('uint%i', TIF.BitsPerSample(1));
    case 2
        classname = sprintf('int%i', TIF.BitsPerSample(1));
    case 3
        if ( TIF.BitsPerSample(1) == 32 )
            classname = 'single';
        else
            classname = 'double';
        end
    otherwise
        error('unsuported TIFF sample format %i', TIF.SampleFormat);
end

image = zeros(height,width,Nchannels,classname);
for c=channels
    % Read the strips
    [image_temp,TIF] = read_strip(width, c, classname,TIF);
    % Extract valid part of data if needed
    if ~all(size(image_temp) == [width height]),
        image_temp = image_temp(1:width, 1:height);
        warning('tiffread2:Crop','Cropping data: found more bytes than needed');
    end
    image(:,:,c) = image_temp;
end


%% ============= sub-functions to read a channel strip ===================

function [strip,TIF] = read_strip(width, stripCnt, classname,TIF)
    StripLength = TIF.StripByteCounts(stripCnt) ./ TIF.BytesPerSample(stripCnt);
    status = fseek(TIF.file, TIF.StripOffsets(stripCnt), 'bof');
    if status == -1
        error('invalid file offset (error on fseek)');
    end
    classnameconv = ([classname '=>' classname]); %read from type to type
    bytes = fread(TIF.file, StripLength, classnameconv, TIF.BOS);
    if any( length(bytes) ~= StripLength )
        error('End of file reached unexpectedly.');
    end
    strip = reshape(bytes,width,StripLength/width);
