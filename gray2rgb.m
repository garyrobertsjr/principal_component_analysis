function [ rgbTlo ] = gray2rgb( Tlo )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Convert greyscale image to RGB
    rgbTlo = repmat(Tlo,[1 1 3]);
    rgbTlo = cat(3,Tlo,Tlo,Tlo);

end

