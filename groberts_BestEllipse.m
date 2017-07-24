%%
% Read Image
image = imread('fruit1.bmp');
color_image=image;
figure;
imshow(image);
title('Original Image');
%%
%Double Threshold
[height, width, channels] = size(image);

if(channels==3)
    image=rgb2gray(image);
end

[Iprime, rgbTlo] = doubleThresh(image);
figure;
imshow(Iprime);
title('Treshold -- Foreground Background separated');

[L,num_components,labels] = ConnectedComponentsByRepeatedFloodfill(Iprime);


u_out = uint8(L);
figure;
imshow(u_out)
A=rgb2gray(u_out);
title('Separate Components');

%% Your algorithm to draw axes under
% Draw the axes for banans only. To find the label of bananas in separate components image, use the figures data cursor option 
t_sf = [];
sf = [];

figure;
imshow(image);
title('Fruits Classified -- Annotated w/ Colors');

% Identify Each Object
for c=1:num_components
    
    % Find pixels belonging to object
    m00=0;m01=0;m10=0;m11=0;m20=0;m02=0;
    for i=1:height
        for j=1:width
            % If pixel belong to object, aggregate moment
            if(L(i,j) == labels(c))
                m01 = m01 + i.^0 * j.^1;
                m00 = m00 + i.^0 * j.^0;
                m10 = m10 + i.^1 * j.^0;
                m11 = m11 + i.^1 * j.^1;
                m20 = m20 + i.^2 * j.^0;
                m02 = m02 + i.^0 * j.^2;
            end
        end
    end
    
    % Centers
    xc = m10/m00;
    yc = m01/m00;
    
    % Catch nan conditions
    if(isnan(xc) || isnan(yc))
        break;
    end
    
    u00 = m00;
    u11 = m11 - yc*m10;
    u20 = m20 - xc*m10;
    u02 = m02 - yc*m01;
    
    C = [u20, u11; u11, u02];
    C = (1/u00).*C;
    
    E = eig(C);
    disp(E);
    theta = .5*atan2(2*u11,u20-u02);
    
    % Eccentricity
    ecc = sqrt((E(2)-E(1))/E(2));

    % major axis
    major_x = [xc - cos(theta)*sqrt(E(2)), xc, xc + cos(theta)*sqrt(E(2))];
    major_y = [yc - sin(theta)*sqrt(E(2)), yc, yc + sin(theta)*sqrt(E(2))]; 
    hold on;

    % minor axis
    minor_x = [xc + cos(-1.5708+theta)*sqrt(E(1)), xc, xc + cos(1.5708+theta)*sqrt(E(1))];
    minor_y = [yc + sin(-1.5708+theta)*sqrt(E(1)), yc, yc + sin(1.5708+theta)*sqrt(E(1))];
    fprintf('Theta %0.5g\n', theta);
    
    % Display center of object
    fprintf('Object: %d\n Intensity: %d\n Center: (%0.5g,%0.5g)\n Eccentricity:%0.5g\n', c, labels(c), xc, yc, ecc);
    disp(m11);
    
    % Classify Bananas
    len = sqrt(power(major_x(3)-major_x(1),2)+power(major_y(3)-major_y(1),2));
    if(ecc >.94 && ecc < .98 && ~isnan(ecc) && len > (width)/5)
        line(major_y,major_x, 'Color', 'y', 'LineWidth', 2);
        line(minor_y,minor_x, 'Color', 'y', 'LineWidth', 2);
    end
    
    % Find spherical fruits
    if(ecc >0 && ecc <.6 && ~isnan(ecc))
        %line(major_y,major_x, 'Color', 'r');
        %line(minor_y,minor_x, 'Color', 'r');
        t_sf = cat(1, t_sf, [major_x, major_y, minor_x, minor_y]);
    end
end

%Compute radius of each spherical fruit
t_radii = [];
f_sf = [];
n_sf = size(t_sf);
radii = [];

for f=1:n_sf(1)
    t_radii = [t_radii sqrt(power(t_sf(f,3)-t_sf(f,1),2)+power(t_sf(f,6)-t_sf(f,4),2))];
end

n_rad = size(t_radii);
r_mean = mean(t_radii);
r_std = std(t_radii);

for f=1:n_sf(1)
    if(t_radii(f) >= r_mean-r_std)
        radii = [radii t_radii(f)];
        %sf = [sf t_sf(f)];
        sf = cat(1, sf, t_sf(f,:));
    end
end
n_rad = size(radii);
n_sf = size(sf);
for f=1:n_sf(1)
    if(radii(f) >= mean(radii)-.01*std(radii))
        line([sf(f,4), sf(f,5), sf(f,6)],[sf(f,1), sf(f,2), sf(f,3)], 'Color', [1,.5,0], 'LineWidth', 2);
        line([sf(f,10), sf(f,11), sf(f,12)],[sf(f,7), sf(f,8), sf(f,9)], 'Color', [1,.5,0], 'LineWidth', 2);
    else
        line([sf(f,4), sf(f,5), sf(f,6)],[sf(f,1), sf(f,2), sf(f,3)], 'Color', 'r', 'LineWidth', 2);
        line([sf(f,10), sf(f,11), sf(f,12)],[sf(f,7), sf(f,8), sf(f,9)], 'Color', 'r', 'LineWidth', 2);
    end
end
