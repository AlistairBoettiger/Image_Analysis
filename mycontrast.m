% Alistair Boettiger

function Iout = mycontrast(I,maxp,minp)

% I is an MxN (2D) image file, can be uint8, uint16, single or double
% maxp is between [0,1], the percent of pixels to saturate
% minp is between [0,1], the precent of pixels to send to 0

s = 1-maxp;
v = sort(I(:));  
l = length(v);

c = class(I);
if strcmp(c,'uint8')==1
    m = 255;
elseif strcmp(c,'uint16')==1
    m = 2^16;
elseif strcmp(c,'single') == 1 || strcmp(c,'double') == 1
    m = max(v);
end

if sum(I(:))~= 0
    o1 = double(v(1+round(minp*l)))/m;
    o2 = double(v(round(s*l)))/m;    
    while o2 == 0
        maxp = maxp*.1;
        s = 1-maxp;
        o2 = double(v(round(s*l)))/m;
    end

    Iout = imadjust(I,[o1,o2],[0,1]); 
else
    Iout = I;
end


%  figure(1); clf; subplot(1,2,1); imagesc(I);
% subplot(1,2,2); imagesc(Iout);