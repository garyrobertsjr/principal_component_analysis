% Author: Alex Murray
% doubleThresh.m: Thresholds an image according to the top 80th percentile
% and lowest 7.5 percentile. This separated items in the foreground and
% background for futher image processing algorithms. Additionally fills
% values higher than the high threshold with white, producing a binary
% image
% Input: Greyscale image
% Output: Thresholded binary image


function [Iprime, rgbTlo] = doubleThresh(img)

    % Get dimensions of original image
    [m,n] = size(img);
    
    % Transform 2-D array into 1-D for sorting purposes
    d1Image = img(:);
    d1Image = d1Image(:);
    
    % Set percentile of thresholding
    hi = 0.80;
    lo = 0.15;

    % Initialize output image
    Iprime = zeros(m,n,3);
    
    % Sort 1-D image
    srt = sort(d1Image);
    
    % Get length of 1-D array
    sz = size(srt);
    
    % Get high and low threshold values 
    tLow = srt(round(lo*sz(1)));
    tHigh = srt(round(hi*sz(1)));
    
    % Get rid of values less than the low threshold
    Tlo=zeros(m,n);
    Thi=zeros(m,n);
    
    for x=1:m
        for y=1:n
            if(img(x,y)>tHigh)
                Thi(x,y)=img(x,y);
            end
            if(img(x,y)>tLow)
                Tlo(x,y)=img(x,y);
            end
        end
    end
    
    
    %Convert greyscale image to RGB
    rgbTlo = repmat(Tlo,[1 1 3]);
    rgbTlo = cat(3,Tlo,Tlo,Tlo);
    rgbThi = repmat(Thi,[1 1 3]);
    rgbThi = cat(3,Thi,Thi,Thi);
    % Set output color to white
    color = [255, 255, 255];
    
    % Go pixel by pixel
    for i=1:m
        for j = 1:n
            % If pixel value is greater than the threshold, then Floodfill 
            % output image at (i,j) with white for an effective threshold
            if(~(rgbThi(i,j,1)==0 && rgbThi(i,j,2)==0 && rgbThi(i,j,3)==0))
                rgbThi=FloodFill(rgbThi,i,j,[0,0,0]);
                Iprime = FloodFillSeparate(rgbTlo,Iprime,i,j,color);
            end
        end
    end
end