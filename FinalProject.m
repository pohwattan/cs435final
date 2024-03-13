% PLEASE CHANGE THIS CODE IT IS UGLY

clear all;
close all;

im = double(imread("first.jpg"))/255;
im2 = double(imread("second.jpg"))/255;

% PREPROCESSING - Rotate images to vertical orientation
im = imrotate(im, -90);
im2 = imrotate(im2, -90);

[im, im2, points, points2] = point_correspondences(im, im2);
im = uint8(im*255);
im2 = uint8(im2*255);
figure(1);
imshow(im);
figure(2);
imshow(im2);


% 1 (10 points) Hard Coding Point Correspondences
% TODO ~ IMPLEMENT DIFFERENT METHOD TO ADDING THE CIRCLES IN THE IMAGES
function [im, im2, points, points2] = point_correspondences(im, im2)
    points = [1390 135; 2177 875; 1205 2045; 2643 2983];
    points2 = [665 86; 1535 947; 473 2032; 1839 2995];
    radius = 20;
    colors = {'red', 'yellow', 'blue', 'green'};
    for i = 1:size(points,1)
        im = insertShape(im, 'FilledCircle', [points(i,:) radius], 'Color', colors{i});
        im2 = insertShape(im2, 'FilledCircle', [points2(i,:) radius], 'Color', colors{i});
    end
end