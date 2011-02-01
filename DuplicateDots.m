

%% Description
% Find replicate dots that are present in the above layer.

% inds_out is a vector of indices at the rescaled resolution 2048/scale
% which contains indices of all the dots in layer j that are not duplicated
% in layer j-1, where duplication is defined as being overlapped within
% 'scale' number of pixels.

function inds_out = DuplicateDots(D1,D2,scale,h,w)


        if isempty(D1)
            [R2,inds2] = vect2rast(D2(:,1),D2(:,2),h,w,scale);
            inds_out = inds2;
        else
 
%         % Plotting for troubleshooting: circles    
%         figure(1); clf; plot(D1(:,1),D1(:,2),'ro'); hold on;
%         plot(D2(:,1),D2(:,2),'go');


       % h = 2048; w = 2048; % localize to single pixel (50 nm) resolution
        % could do Gaussian fitting and localize to better than this).  
        [R1,inds1] = vect2rast(D1(:,1),D1(:,2),h,w,scale);
        [R2,inds2] = vect2rast(D2(:,1),D2(:,2),h,w,scale);


%        % Labeling and plotting overlapped pixels:
%        figure(2); clf; imagesc(R1+2*R2);
        overlap = R1+2*R2; 

% % Plotting just Overlapped pixels
%         [h1,w1] = size(overlap);
%         doubles = uint8(zeros(h1,w1)); 
%         doubles(overlap==3) = 1; % map of doubles;
%         figure(2); clf; imshow(255*doubles);

    % remove overlapped pixels from pixel id list.  
        inddoubles = find(overlap==3);  % linear indices of doubles 
        inds_out = setdiff(inds2,inddoubles);

        end

