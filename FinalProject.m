clear variables;
close all;

im = imread("first.jpg");
im2 = imread("second.jpg");

% PREPROCESSING - Rotate images to vertical orientation
im = imrotate(im, -90);
im2 = imrotate(im2, -90);

%%%%%%%%%%%%%%
% QUESTION 1 %
%%%%%%%%%%%%%%
[imMarked, im2Marked, points, points2] = point_correspondences(im, im2);
figure('Name','Point Correspondence Side-by-Side', 'FileName','PointCorrespondence.jpg');
imshowpair(imMarked,im2Marked,'montage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 2 by Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transformationMatrix = FindTransformationMatrixWithPoints(points, points2);
stitchedImage = StitchImages(im, im2, transformationMatrix);
figure('Name','Stitched Image', 'FileName','StitchedImage.jpg');
imshow(stitchedImage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 3 by Nick Pohwat %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pyramids = image_pyramids(im, im2);
for i = 1:length(pyramids)
    figure('Name', ['Scale-Space Image Pyramid ' num2str(i)], 'FileName', ['ImagePyramid' num2str(i) '.jpg']);
    show_pyramids(pyramids{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 4 by Nick Pohwat and Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[all_extrema_image, pruned_extrema_image, remaining_extrema] = local_maximas(pyramids{1});
figure('Name','Finding the Local Maximas 1', 'FileName','LocalMaximas1.jpg');
imshowpair(all_extrema_image,pruned_extrema_image,'montage');
[all_extrema_image2, pruned_extrema_image2, remaining_extrema2] = local_maximas(pyramids{2});
figure('Name','Finding the Local Maximas 2', 'FileName','LocalMaximas2.jpg');
imshowpair(all_extrema_image2,pruned_extrema_image2,'montage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 5 by Nick Pohwat %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [left_images, right_images] = keypoint_matching(im, im2, remaining_extrema, remaining_extrema2);
% figure('Name','Keypoint Description and Matching', 'FileName','KeypointMatching.jpg');
% imshowpair(left_images,right_images,'montage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 6 by Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1 (10 points) Hard Coding Point Correspondences
function [im, im2, points, points2] = point_correspondences(im, im2)
    points = [545 238; 299 510; 665 896; 281 896];
    points2 = [383 258; 115 506; 443 896; 65 918];
    radius = 8;
    colors = {'red', 'yellow', 'blue', 'green'};
    for i = 1:size(points,1)
        im = insertShape(im, 'FilledCircle', [points(i,:) radius], 'Color', colors{i});
        im2 = insertShape(im2, 'FilledCircle', [points2(i,:) radius], 'Color', colors{i});
    end
end

% Get the transformation matrix using points
function homogenousMatrix = FindTransformationMatrixWithPoints(basePoints, primePoints)
    %Note: Points are still in (x,y) notation
    
    homogenousMatrix = zeros(3,3);
    A = [];
    B = [];
    rowNum = 1;
    for pointNum = 1:4
        A(rowNum, :) = [basePoints(pointNum, 1) basePoints(pointNum, 2) 1 0 0 0 0 0 0];
        B(rowNum, :) = [0 0 0 0 0 0 (basePoints(pointNum, 1) * primePoints(pointNum, 1)) (basePoints(pointNum, 2) * primePoints(pointNum, 1)) primePoints(pointNum, 1)];
        rowNum = rowNum + 1;
        A(rowNum, :) = [0 0 0 basePoints(pointNum, 1) basePoints(pointNum, 2) 1 0 0 0];
        B(rowNum, :) = [0 0 0 0 0 0 (basePoints(pointNum, 1) * primePoints(pointNum, 2)) (basePoints(pointNum, 2) * primePoints(pointNum, 2)) primePoints(pointNum, 2)];
        rowNum = rowNum + 1;
    end
    At = A.';
    Bt = B.';
    [~, ~, V] = svd((At * A) + (Bt * B) - (Bt * A) - (At * B));
    m = V(:,end);
    count = 1;
    for y = 1:3
        for x = 1:3
            homogenousMatrix(y,x) = m(count,:);
            count = count + 1;
        end
    end
    
    %For sanity checking, transform a random point and see if it is correct!
    testPointNum = randi(4);
    testPoint = basePoints(testPointNum,:);
    testPointPrime = primePoints(testPointNum, :);
    transformTestPoint = TransformPoint(testPoint, homogenousMatrix);
    isTransformationCorrect = isequal(testPointPrime, transformTestPoint);
    if (~isTransformationCorrect) 
        disp("Something went wrong finding the transformation matrix!");
    end
end

function transformedPoints = TransformMultiplePoints(points, transformationMatrix)
    transformedPoints = zeros(size(points, 1), size(points, 2));
    for pointNum = 1:size(points,1)
        transformedPoints(pointNum, :) = TransformPoint(points(pointNum,:), transformationMatrix);
    end
end

function transformedPoint = TransformPoint(point, transformationMatrix)
    %Assuming points are coming in as (x,y)
    homogenousPoint = [
            point(1);
            point(2);
            1
        ];
    homogenousTransformedPoint = transformationMatrix * homogenousPoint;
    transformedPoint = [round(homogenousTransformedPoint(1,:)/homogenousTransformedPoint(3,:)), round(homogenousTransformedPoint(2,:)/homogenousTransformedPoint(3,:))];
end

function stitchedImage = StitchImages(baseImage, transformedImage, transformationMatrix)
    %Points are still in (x,y) format
    transformedImageCornerPoints = [
        1 1;                                                %Top Left
        1 size(transformedImage, 1);                        %Bottom Left
        size(transformedImage, 2) 1;                        %Top Right
        size(transformedImage, 2) size(transformedImage, 1) %Bottom Right
    ];
    transformedImageTransformedCornerPoints = TransformMultiplePoints(transformedImageCornerPoints, inv(transformationMatrix));
    minTransformedY = min(transformedImageTransformedCornerPoints(1,2), transformedImageTransformedCornerPoints(3,2));
    maxTransformedY = max(transformedImageTransformedCornerPoints(2,2), transformedImageTransformedCornerPoints(4,2));
    minTransformedX = min(transformedImageTransformedCornerPoints(1,1), transformedImageTransformedCornerPoints(2,1));
    maxTransformedX = max(transformedImageTransformedCornerPoints(3,1), transformedImageTransformedCornerPoints(4,1));
    % transformedImageTransformedHeight = maxTransformedY - minTransformedY;
    % transformedImageTransformedWidth = maxTransformedX - minTransformedX;
    
    minBaseImageY = 0;
    maxBaseImageY = size(baseImage, 1);
    minBaseImageX = 0;
    maxBaseImageX = size(baseImage, 2);
    % baseImageHeight = size(baseImage, 1);
    baseImageWidth = size(baseImage, 2);
    
    %Setup canvas for the stitched image
    minStitchedImageY = min(minBaseImageY, minTransformedY);
    maxStitchedImageY = max(maxBaseImageY, maxTransformedY);
    minStitchedImageX = min(minBaseImageX, minTransformedX);
    maxStitchedImageX = max(maxBaseImageX, maxTransformedX);
    stitchedImageHeight = maxStitchedImageY - minStitchedImageY;
    stitchedImageWidth = maxStitchedImageX - minStitchedImageX;
    stitchedImage = uint8(zeros(stitchedImageHeight, stitchedImageWidth, 3));

    for y = 1:stitchedImageHeight
        clc;
        disp("Stitching Progress = " + round((y/stitchedImageHeight)*100) + "%");
        for x = 1:stitchedImageWidth
            if (minTransformedX < 1)
                baseImageX = x + minTransformedX;
            else
                baseImageX = x;
            end

            if (minTransformedY < 1)
                baseImageY = y + minTransformedY;
            else
                baseImageY = y;
            end

            baseImagePoint = [baseImageX, baseImageY];
            transformedImagePoint = TransformPoint(baseImagePoint, transformationMatrix);
            transformedImageY = transformedImagePoint(2);
            transformedImageX = transformedImagePoint(1);
            
            pointExistsInBaseImage = PointInImage(baseImagePoint, baseImage);
            pointExistsInTransformedImage = PointInImage(transformedImagePoint, transformedImage);

            if (~pointExistsInBaseImage)
                baseImagePointColor = [0, 0, 0];
            else
                baseImagePointColor = [baseImage(baseImageY, baseImageX, 1) baseImage(baseImageY, baseImageX, 2) baseImage(baseImageY, baseImageX, 3)];
            end

            if (~pointExistsInTransformedImage)
                transformedImagePointColor = [0, 0, 0];
            else
                transformedImagePointColor = [transformedImage(transformedImageY, transformedImageX, 1) transformedImage(transformedImageY, transformedImageX, 2) transformedImage(transformedImageY, transformedImageX, 3)];
            end
            
            if (pointExistsInBaseImage && pointExistsInTransformedImage)
                stitchedImage(y,x,:) = uint8(BlendPointColors(baseImagePointColor, transformedImagePointColor, (1 - baseImageX/baseImageWidth)));
            elseif (pointExistsInTransformedImage)
                stitchedImage(y,x,:) = uint8(transformedImagePointColor);
            else
                stitchedImage(y,x,:) = uint8(baseImagePointColor);
            end
        end
    end
end 

function existsInImage = PointInImage (point, image)
    existsInImage = (point(1) > 0) && (point(1) <= size(image,2)) && (point(2) > 0) && (point(2) <= size(image,1));
end

function blendedPointColors = BlendPointColors (point1Colors, point2Colors, alpha)
    blendedPointColors = [point1Colors(1)*alpha+point2Colors(1)*(1-alpha) point1Colors(2)*alpha+point2Colors(2)*(1-alpha) point1Colors(3)*alpha+point2Colors(3)*(1-alpha)];
end

% 3 (10 points) Create Scale-Space Image Pyramids
function pyramids = image_pyramids(varargin)
    filterSize = @(sigma) ceil(3*sigma)*2+1;
    m_scales = 5;
    n_octaves = 4;

    pyramids = cell(length(varargin), 1);
    for i = 1:length(varargin)
        im = varargin{i};
        if size(im, 3) == 3
            im = (0.2989*im(:,:,1)) + (0.5870*im(:,:,2)) + (0.1140*im(:,:,3));
        end

        pyramid = cell(n_octaves, m_scales);
        for n = 1:n_octaves
            k = sqrt(2);
            for m = 1:m_scales
                sigma = 2^(n-1) * k^(m-1) * 1.6;
                pyramid{n, m} = imgaussfilt(im, sigma, 'FilterSize', filterSize(sigma));
            end
            im = im(1:2:end, 1:2:end);
        end
        pyramids{i} = pyramid;
    end
end

function show_pyramids(varargin)
    for i = 1:length(varargin)
        pyramid = varargin{i};
        [height, width] = size(pyramid);
        for n = 1:height
            for m = 1:width
                subplot(height, width, m + ((n-1) * width));
                imshow(pyramid{n, m});
                axis on;
            end
        end
    end
end

% 4 (10 points) Finding the Local Maximas
function [all_extrema_image, pruned_extrema_image, pruned_extrema_bits] = local_maximas(pyramid)
    all_extrema_image = pyramid{1, 1};
    pruned_extrema_image = pyramid{1, 1};
    location_multiplier = @(octave) 2^(octave - 1);
    all_extremas_bits = zeros(size(all_extrema_image));

    for octave = 1:size(pyramid, 1)
        for scale = 2:size(pyramid, 2) - 1
            scale_above = pyramid{octave,scale-1};
            scale_current = pyramid{octave,scale};
            scale_below = pyramid{octave,scale+1};
            for height = 2:size(scale_current,1)-1
                for width = 2:size(scale_current,2)-1
                    window = [scale_above(height-1:height+1, width-1:width+1);
                              scale_current(height-1:height+1, width-1:width+1);
                              scale_below(height-1:height+1, width-1:width+1)];

                    max_val = max(window(:));
                    min_val = min(window(:));

                    if scale_current(height, width) == max_val || scale_current(height, width) == min_val
                        scaled_height = height * location_multiplier(octave);
                        scaled_width = width * location_multiplier(octave);
                        all_extremas_bits(scaled_height, scaled_width) = 1;
                    end
                end
            end
        end
    end
    [extrema_height, extrema_width] = find(all_extremas_bits == 1);
    all_extrema_image = insertShape(all_extrema_image, 'Circle', [extrema_width, extrema_height, ones(length(extrema_width), 1)*5], 'Color', 'red');

    pruned_extrema_bits = all_extremas_bits;

    [extrema_height, extrema_width] = find(pruned_extrema_bits == 1);
    edges = edge(pruned_extrema_image, 'Canny');
    for i = 1:length(extrema_height)
        if edges(extrema_height(i), extrema_width(i))
            pruned_extrema_bits(extrema_height(i), extrema_width(i)) = 0;
        end
    end
    
    border_distance = 9;
    pruned_extrema_bits(1:border_distance, :) = 0;
    pruned_extrema_bits(end-border_distance+1:end, :) = 0;
    pruned_extrema_bits(:, 1:border_distance) = 0;
    pruned_extrema_bits(:, end-border_distance+1:end) = 0;

    std_dev_threshold = 0.7;
    patch_size = 9;
    patch_radius = floor(patch_size / 2);
    [extrema_height, extrema_width] = find(pruned_extrema_bits == 1);
    for i = 1:length(extrema_height)
        patch = pruned_extrema_image(extrema_height(i) - patch_radius:extrema_height(i) + patch_radius,extrema_width(i) - patch_radius:extrema_width(i) + patch_radius);
        patch = double(patch(:));
        if std(patch(:)) < std_dev_threshold
            pruned_extrema_bits(extrema_height(i), extrema_width(i)) = 0;
        end
    end
    
    [extrema_height, extrema_width] = find(pruned_extrema_bits == 1);
    pruned_extrema_image = insertShape(pruned_extrema_image, 'Circle', [extrema_width, extrema_height, ones(length(extrema_width), 1)*5], 'Color', 'red');
end
