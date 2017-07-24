% Author: Alexander Murray
% ConnectedComponentsByRepeatedFloodfill.m: Used for finding connected
% components in an image. Used for differentiating between objects in an image
% Input: RGB Image (Works better when thresholded)
% Output: Labeled Image with a greyscale pixel value labels
% Requires: Floodfill.m, pop.m, push.m

function [L,num_components,labels] = ConnectedComponentsByRepeatedFloodfill(image)

    % Get dimensions of input image
    [h,w,c] = size(image);
    % Initializing output image
    L = double(zeros(h,w,c));
    num_components=0;
    % Initialize label to black
    next_label = [0 0 0];
    
    % Go pixel by pixel
    for i = 1:h
        for j = 1:w
            % If pixel is empty, then perform operations
            if(L(i,j) == 0)
                
                % Floodfill image at (i,j) with incremented label
                L = FloodFillSeparate(image, L, i, j, next_label);
                num_components=num_components+1; 
                labels(num_components)=next_label(1) + 10;
                % Increment Label
                next_label(1) = next_label(1) + 10;
                next_label(2) = next_label(2) + 10;
                next_label(3) = next_label(3) + 10;
            end
        end
    end


end