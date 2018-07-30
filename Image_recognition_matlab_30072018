clc;
clear;
clear all;

%%image resmi okuyor

image_original = imread('image1.jpeg');

%% image gray yapar

image_gray = rgb2gray(image_original);

%% image sobel yapar

image_sobel = edge(image_gray , 'sobel');

%% Binary image
BW = im2bw(image_sobel,0.6);
% BW = ~BW;

%% boşlukları doldurur
% % se=strel('disk',10);
% % bw=imclose(image_BW,se);


se = strel('disk' , 10) %disk biçiminde yapısal element oluşturuyoruz
bw = imclose(BW,se) %iç kısmındaki boşluklar kayboldu
BW = imfill(bw,'holes');

%% tekrar gray yapar
image_BW = im2double(BW)

%% gausfilter method
image_gaussfilter = imgaussfilt(image_BW,0.75);
image_gaussfilter = imgaussfilt(image_gaussfilter,1.65);

% % cornermetric method
% image_corner = cornermetric(image_just)

%%
image_BW = im2double(image_gaussfilter)


%% imadjust method
image_just = imadjust(image_BW);

%% konkav yapar
image_conkav = bwconvhull(image_just,'objects');
image_conkav = im2double(image_conkav)


%% show
figure
imshow(image_original)
figure
imshow(image_just);

%% image small yapar
image_small = imresize(image_just,[400 600]);


% %%imgradientxy methodu w(x,y)
% [Gx,Gy] = imgradientxy(image_small);
% figure
% imshowpair(Gx,Gy, 'montage');


%% function yap burada

% %% compute x and y derivated  image
% Ix = Gx * image_small;
% Iy = Gy * image_small;
% 
% %% comute the sum of product of derivated pixel
% Ix2 = Ix * Ix;
% Iy2 = Iy * Iy;
% Ixy = Ix * Iy; 
% 
% Sx2 = [Gx,Gy] * Ix2;



k = 0.06;
Threshold = 100000;
sigma = 2;
halfwid = sigma * 3;

[xx, yy] = meshgrid(-halfwid:halfwid, -halfwid:halfwid);

Gxy = exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));

Gx = xx .* exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));
Gy = yy .* exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));

numOfRows = size(image_small, 1);
numOfColumns = size(image_small, 2);

% 1) Compute x and y derivatives of image
Ix = conv2(Gx, image_small);
Iy = conv2(Gy, image_small);

size(Ix);

% 2) Compute products of derivatives at every pixel
Ix2 = Ix .^ 2;
Iy2 = Iy .^ 2;
Ixy = Ix .* Iy;

% 3)Compute the sums of the products of derivatives at each pixel
Sx2 = conv2(Gxy, Ix2);
Sy2 = conv2(Gxy, Iy2);
Sxy = conv2(Gxy, Ixy);

im = zeros(numOfRows, numOfColumns);
for x = 1:numOfRows,
   for y = 1:numOfColumns,
       %x,y
       % 4) Define at each pixel(x, y) the matrix H
       H = [Sx2(x, y) Sxy(x, y); Sxy(x, y) Sy2(x, y)];
       
       % 5) Compute the response of the detector at each pixel
       R = det(H) - k * (trace(H) ^ 2);
       
       % 6) Threshold on value of R
       if (R > Threshold)
          im(x, y) = R;
       end
   end
end

% 7) Compute nonmax suppression
output = im > imdilate(im , [1 1 1; 1 0 1; 1 1 1]);
output_1 = output + image_small;
% figure, imshow(image_small);
figure
imshow(image_small);
figure
imshow(output)





%%
% figure
% corners = detectHarrisFeatures(image_original);
% imshow(image_original)
% hold on;
% plot(corners.selectStrongest(50));
% title('detectHarrisFeatures')

%% köşe bulma
 
% C = corner(image_conkav,1);  %// specifying maximum no. of corners
% 
% figure;
% imshow(image_conkav);
% hold on
% scatter(C(:,1),C(:,2),50,'filled');



