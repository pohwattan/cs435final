clear variables;
close all;

im = imread("first.jpg");
im2 = imread("second.jpg");

% PREPROCESSING - Rotate images to vertical orientation
im = imrotate(im, -90);
im2 = imrotate(im2, -90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 1 by Nick Pohwat and Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

<<<<<<< HEAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 3 by Nick Pohwat %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 3 by Nick Pohwat and Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>>>>>>> 90729ca (start part 5)
[pyramids, DoGs] = image_pyramids(im, im2);
for i = 1:length(pyramids)
    figure('Name', ['Scale-Space Image Pyramid ' num2str(i)], 'FileName', ['ImagePyramid' num2str(i) '.jpg']);
    show_pyramids(pyramids{i});
    figure('Name', ['DoG Pyramid ' num2str(i)], 'FileName', ['DoGPyramid' num2str(i) '.jpg']);
    show_pyramids(DoGs{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 4 by Nick Pohwat and Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
<<<<<<< HEAD
[all_extrema_image, pruned_extrema_image, remaining_extrema] = local_maximas(DoGs{1}, pyramids{1}{1, 1});
figure('Name','Finding the Local Maximas 1', 'FileName','LocalMaximas1.jpg');
imshowpair(all_extrema_image,pruned_extrema_image,'montage');
[all_extrema_image2, pruned_extrema_image2, remaining_extrema2] = local_maximas(DoGs{2}, pyramids{2}{1, 1});
=======
[all_extrema_image, pruned_extrema_image, remaining_extrema] = local_maximas(DoGs{1});
figure('Name','Finding the Local Maximas 1', 'FileName','LocalMaximas1.jpg');
imshowpair(all_extrema_image,pruned_extrema_image,'montage');
[all_extrema_image2, pruned_extrema_image2, remaining_extrema2] = local_maximas(DoGs{2});
>>>>>>> 90729ca (start part 5)
figure('Name','Finding the Local Maximas 2', 'FileName','LocalMaximas2.jpg');
imshowpair(all_extrema_image2,pruned_extrema_image2,'montage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 5 by Nick Pohwat and Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[left_images, right_images] = keypoint_matching(im, im2, remaining_extrema, remaining_extrema2);
figure('Name','Keypoint Description and Matching', 'FileName','KeypointMatching.jpg');
imshowpair(left_images,right_images,'montage');

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
function [pyramids, DoGs] = image_pyramids(varargin)
    kernel = @(sigma) ceil(3*sigma)*2+1;
    m_scales = 5;
    n_octaves = 4;

    pyramids = cell(length(varargin), 1);
    DoGs = cell(length(varargin), 1);
    for i = 1:length(varargin)
        im = double(varargin{i}) / 255;
        if size(im, 3) == 3
            im = (0.2989*im(:,:,1)) + (0.5870*im(:,:,2)) + (0.1140*im(:,:,3));
        end

        pyramid = cell(n_octaves, m_scales);
        DoG = cell(n_octaves, m_scales-1);
        for n = 1:n_octaves
            for m = 1:m_scales
                sigma = 2^(n-1) * sqrt(2)^(m-1) * 1.6;
                pyramid{n, m} = imgaussfilt(im, sigma, 'FilterSize', kernel(sigma));
                if m >= 2
                    DoG{n,m-1} = pyramid{n, m-1} - pyramid{n, m};
                end
            end
            im = im(1:2:end, 1:2:end);
        end
        pyramids{i} = pyramid;
        DoGs{i} = DoG;
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

function [extremaValue, extremaLocation] = FindMatrixExtrema(inputMatrix, mode)
    matrixHeight = size(inputMatrix, 1);
    matrixWidth = size(inputMatrix, 2);
    matrixDepth = size(inputMatrix, 3);
    extremaValue = inputMatrix(2,2,2);
    extremaLocation = [2,2,2];
    for y = 1:matrixHeight
        for x = 1:matrixWidth
            for z = 1:matrixDepth
                if (mode == "max")
                    if inputMatrix(x,y,z) > extremaValue
                        extremaValue = inputMatrix(x,y,z);
                        extremaLocation = [x,y,z];
                    end
                elseif (mode == "min")
                    if inputMatrix(x,y,z) < extremaValue
                        extremaValue = inputMatrix(x,y,z);
                        extremaLocation = [x,y,z];
                    end
                else
                    disp("Something went wrong with MinMaxExtrema, could not recognize mode!");
                end
            end
        end
    end
end

function isCurrentPointExtrema = IsCurrentPointExtrema(inputMatrix, mode)
    matrixHeight = size(inputMatrix, 1);
    matrixWidth = size(inputMatrix, 2);
    matrixDepth = size(inputMatrix, 3);
    currentValue = inputMatrix(2,2,2);
    isCurrentPointExtrema = false;
    for y = 1:matrixHeight
        for x = 1:matrixWidth
            for z = 1:matrixDepth
                if (x == 2 && y == 2 && z == 2)
                    continue
                end

                if (mode == "max")
                    if inputMatrix(x,y,z) > currentValue
                        isCurrentPointExtrema = false;
                        return;
                    end
                elseif (mode == "min")
                    if inputMatrix(x,y,z) < currentValue
                        isCurrentPointExtrema = false;
                        return;
                    end
                else
                    disp("Something went wrong with MinMaxExtrema, could not recognize mode!");
                end
            end
        end
    end
    isCurrentPointExtrema = true;
end

% 4 (10 points) Finding the Local Maximas
function [all_extrema_image, pruned_extrema_image, pruned_extrema_bits] = local_maximas(DoGPyramid, baseImage)
    all_extrema_image = DoGPyramid{1, 1};
    pruned_extrema_image = DoGPyramid{1, 1};
    location_multiplier = @(octave) 2^(octave - 1);
    all_extremas_bits = zeros(size(all_extrema_image));

    for octave = 1:size(DoGPyramid, 1)
        for scale = 2:size(DoGPyramid, 2) - 1
            scale_above = DoGPyramid{octave,scale-1};
            scale_current = DoGPyramid{octave,scale};
            scale_below = DoGPyramid{octave,scale+1};
            for height = 2:size(scale_current,1)-1
                for width = 2:size(scale_current,2)-1
                    window = zeros(3,3,3);
                    window(:,:,1) = scale_above(height-1:height+1, width-1:width+1);
                    window(:,:,2) = scale_current(height-1:height+1, width-1:width+1);
                    window(:,:,3) = scale_below(height-1:height+1, width-1:width+1);
                    
                    
                    isCurrentPointMax = IsCurrentPointExtrema(window, "max");
                    isCurrentPointMin = IsCurrentPointExtrema(window, "min");

                    % imageSpace_max_loc_x = (max_loc(2) - 2) + width * location_multiplier(octave);
                    % imageSpace_max_loc_y = (max_loc(1) - 2) + height * location_multiplier(octave);
                    % all_extremas_bits(imageSpace_max_loc_y, imageSpace_max_loc_x) = 1;
                    % 
                    % imageSpace_min_loc_x = (min_loc(2) - 2) + width * location_multiplier(octave);
                    % imageSpace_min_loc_y = (min_loc(1) - 2) + height * location_multiplier(octave);
                    % all_extremas_bits(imageSpace_min_loc_y, imageSpace_min_loc_x) = 1;
                    if (isCurrentPointMax || isCurrentPointMin)
                        scaled_height = height * location_multiplier(octave);
                        scaled_width = width * location_multiplier(octave);
                        all_extremas_bits(scaled_height, scaled_width) = 1;
                    end
                end
            end
        end
    end
    [extrema_height, extrema_width] = find(all_extremas_bits == 1);
    all_extrema_image = insertShape(baseImage, 'Circle', [extrema_width, extrema_height, ones(length(extrema_width), 1)*5], 'Color', 'red');

    pruned_extrema_bits = all_extremas_bits;

    edgeBits = edge(pruned_extrema_image, 'Canny');
    for i = 1:length(extrema_height)
        if (edgeBits(extrema_height(i), extrema_width(i)) == 1)
            pruned_extrema_bits(extrema_height(i), extrema_width(i)) = 0;
        end
    end
    
    border_distance = 9;
    [extrema_height, extrema_width] = find(pruned_extrema_bits == 1);
    for i = 1:length(extrema_height)
        if extrema_height(i) <= border_distance ||...
            extrema_width(i) <= border_distance ||...
            extrema_height(i) > size(pruned_extrema_image, 1) - border_distance ||...
            extrema_width(i) > size(pruned_extrema_image, 2) - border_distance
            pruned_extrema_bits(extrema_height(i), extrema_width(i)) = 0;
        end
    end
    % pruned_extrema_bits(1:border_distance, :) = 0;
    % pruned_extrema_bits(end-border_distance+1:end, :) = 0;
    % pruned_extrema_bits(:, 1:border_distance) = 0;
    % pruned_extrema_bits(:, end-border_distance+1:end) = 0;

    std_dev_threshold = 0.03;
    patch_size = 9;
    patch_radius = floor(patch_size / 2);
<<<<<<< HEAD

    [pruned_extrema_height, pruned_extrema_width] = find(pruned_extrema_bits == 1);
    for i = 1:length(pruned_extrema_height)
        patch = pruned_extrema_image(pruned_extrema_height(i) - patch_radius:pruned_extrema_height(i) + patch_radius,pruned_extrema_width(i) - patch_radius:pruned_extrema_width(i) + patch_radius);
=======
    [extrema_height, extrema_width] = find(pruned_extrema_bits == 1);
    for i = 1:length(extrema_height)
        patch = pruned_extrema_image(extrema_height(i) - patch_radius:extrema_height(i) + patch_radius,...
                                     extrema_width(i) - patch_radius:extrema_width(i) + patch_radius);
>>>>>>> 90729ca (start part 5)
        patch = double(patch(:));
        if std(patch(:)) < std_dev_threshold
            pruned_extrema_bits(pruned_extrema_height(i), pruned_extrema_width(i)) = 0;
        end
    end
    
    [pruned_extrema_height, pruned_extrema_width] = find(pruned_extrema_bits == 1);
    pruned_extrema_image = insertShape(baseImage, 'Circle', [pruned_extrema_width, pruned_extrema_height, ones(length(pruned_extrema_width), 1)*5], 'Color', 'red');
end

% 5 (10 points) Keypoint Description and Matching
function [left, right] = keypoint_matching(im, im2, remaining_extrema1, remaining_extrema2)
    descriptors1 = extract_descriptors(im, remaining_extrema1);
    descriptors2 = extract_descriptors(im2, remaining_extrema2);
    left = 0; % Placeholder
    right = 0; % Placeholder
end

function descriptors = extract_descriptors(im, remaining_extrema)
    [keypoint_height, keypoint_width] = find(remaining_extrema == 1);
    num_keypoints = length(keypoint_rows);
    descriptors = zeros(num_keypoints, 243);
    for i = 1:num_keypoints
        patch = im(keypoint_rows(i)-4:keypoint_rows(i)+4, keypoint_cols(i)-4:keypoint_cols(i)+4, :);
        descriptors(i,:) = reshape(patch, [1, 243]);
    end
end
