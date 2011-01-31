
%%                              vect2rast.m
% Alistair Boettiger                                Date Begun: 01/30/11
%
%% Description
% Convert vector data to raster data 
% Returns linear indices of points and Rastermap.  

    % x -- x coordinates of position vector
    % y -- y coordinates of position vector
    % w -- width of raster map
    % h -- height of raster map

%    % Sample data for troubleshooting 
    % x = D2(:,1); y = D2(:,2); 
    % h = 2048; w=h; h*w
    
    function [rasterMap, inds] = vect2rast(x,y,w,h,scale)
    
            x = x/scale; 
            y = y/scale; 
            h = round(h/scale);
            w = round(w/scale); 
            inds = floor(y)+floor(x)*h;  % convert x-y indexing to linear indexing   
            rasterMap = false(h,w); 
            rasterMap(inds) = 1; % raster map of RNA1 positions
            
            % figure;             imshow(rasterMap); 