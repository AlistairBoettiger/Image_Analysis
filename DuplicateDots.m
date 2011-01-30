

%% Description
% Find replicate dots that are present in the above layer.



function [inddoubles] = DuplicateDots(DotData,Zs)
    for j=4:Zs-1
        D1 = DotData{j};
        D2 = DotData{j+1};


        figure(1); clf; plot(D1(:,1),D1(:,2),'ro'); hold on;
        plot(D2(:,1),D2(:,2),'go');


        h = 2055; w = 2055; % localize to single pixel (50 nm) resolution
        % could do Gaussian fitting and localize to better than this).  

        [R1,inds1] = vect2rast(D1(:,1),D1(:,2),h,w,3);
        [R2,inds2] = vect2rast(D2(:,1),D2(:,2),h,w,3);

        figure(2); clf; imagesc(R1+2*R2);
        overlap = R1+2*R2; 

        [h1,w1] = size(overlap);
        doubles = uint8(zeros(h1,w1)); 
        doubles(overlap==3) = 1; % map of doubles;
        figure(2); clf; imshow(255*doubles);

        inddoubles = find(overlap==3);  % linear indices of doubles
    end


% count dots per cell 
    
