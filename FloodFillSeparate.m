% Author: Alex Murray
% FloodFill: Determines area connected to a given node and fills with a
% user specified color. It is used in Microsoft Paint as a "buvket", as 
% well as in Minesweeper to determine cleared pieces (src: Wikipedia)
% Input: RGB Image Matrix, current output image, Pixel coordinate i, j,
% desired color (In RGB format)
% Output: RGB Image Matrix with floodfilled portion

function out = FloodFillSeparate(image, out, i, j, new_color)
    % Get RGB value at (i,j) in the original image
    org_color = [image(i,j,1), image(i,j,2), image(i,j,3)];
    % Get dimensions of input image
    [h,w,c] = size(image);
    
    % Initialize stack for x and y coordinates
    frontierx = zeros(1,h*w);
    frontiery = zeros(1,h*w);
    
    % Initialize stack pointers to zero
    fpx = 0;
    fpy = 0;
   
    % Push pixel coordinates onto stack
    [frontierx, fpx] = push(frontierx,i,fpx);
    [frontiery, fpy] = push(frontiery,j,fpy);
    
    % Give pixel at (i,j) of output image the desired color
    out(i,j,1) = new_color(1);
    out(i,j,2) = new_color(2);
    out(i,j,3) = new_color(3);
   
    % While the stack is not empty
    while (fpx ~= 0 && fpy ~= 0)
        
        % Pop pixel coordinates off stack and adjust 
        [fpx, i] = pop(frontierx, fpx);
        [fpy, j] = pop(frontiery, fpy);
        
        % For 8 surrounding pixels, perform operations
        for k = -1:1
            for l = -1:1
                % Only use values that are within the boundaries of the image
                if(((i+k) <= h) && ((i+k) >= 1) && ((j+l) <= w) && ((j+l) >= 1))
                    % Place color values of the original image and utput
                    % image into RGB format
                    im_color = [image(i+k, j+l, 1), image(i+k, j+l, 2), image(i+k, j+l, 3)];
                    o_color = [out(i+k, j+l, 1), out(i+k, j+l, 2), out(i+k, j+l, 3)];
                    
                    % Use isequal because we are using RGB values
                    if(isequal(im_color,org_color) && ~isequal(new_color, o_color))
                        % Push i+k,j+l onto stack
                        [frontierx, fpx] = push(frontierx, i+k, fpx);
                        [frontiery, fpy] = push(frontiery, j+l, fpy);
                        % Set color of pixels at (i+k, j+l) of output image
                        out(i+k,j+l,1) = new_color(1);
                        out(i+k,j+l,2) = new_color(2);
                        out(i+k,j+l,3) = new_color(3);
                    end
                end
            end
        end
    end
end