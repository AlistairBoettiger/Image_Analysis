
%%                              vect2rast.m
% Alistair Boettiger                                Date Begun: 01/30/11
%
% Convert vector data to raster data 
% 

    function rasterMap = vect2rast(x,y,w,h)
            inds = round(y)+round(x)*h;  % convert x-y indexing to linear indexing   
            rasterMap = false(h,w); 
            rasterMap(inds) = 1; % raster map of RNA1 positions