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
points = [545 238; 299 510; 665 896; 281 896];
points2 = [383 258; 115 506; 443 896; 65 918];
[imMarked, im2Marked, points, points2] = point_correspondences(im, im2, points, points2);
figure('Name','Point Correspondence Side-by-Side', 'FileName','PointCorrespondence.jpg');
imshowpair(imMarked,im2Marked,'montage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 2 by Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transformationMatrix = FindTransformationMatrixWithPoints(points, points2);
stitchedImage = StitchImages(im, im2, transformationMatrix);
figure('Name','Stitched Image', 'FileName','StitchedImage.jpg');
imshow(stitchedImage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 3 by Nick Pohwat and Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pyramids, DoGs] = image_pyramids(im, im2);
for i = 1:length(pyramids)
    figure('Name', ['Scale-Space Image Pyramid ' num2str(i)], 'FileName', ['ImagePyramid' num2str(i) '.jpg']);
    show_pyramids(pyramids{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 4 by Nick Pohwat and Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[all_extrema_image, pruned_extrema_image, remaining_extrema] = local_maximas(DoGs{1}, pyramids{1}{1, 1});
figure('Name','Finding the Local Maximas 1', 'FileName','LocalMaximas1.jpg');
imshowpair(all_extrema_image,pruned_extrema_image,'montage');
[all_extrema_image2, pruned_extrema_image2, remaining_extrema2] = local_maximas(DoGs{2}, pyramids{2}{1, 1});
figure('Name','Finding the Local Maximas 2', 'FileName','LocalMaximas2.jpg');
imshowpair(all_extrema_image2,pruned_extrema_image2,'montage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 5 by Nick Pohwat and Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[canvas, keypoints1, keypoints2, C_union] = keypoint_matching(im, im2, remaining_extrema, remaining_extrema2);
figure('Name','Keypoint Description and Matching', 'FileName','KeypointMatching.jpg');
imshow(canvas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 6 by Nick Pohwat %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[best_canvas, best_transformation_matrix] = auto_stitch(im, im2, keypoints1, keypoints2, C_union);
figure('Name','Best Keypoint Description and Matching', 'FileName','BestKeypointMatching.jpg');
imshow(best_canvas);
bestStitchedImage = StitchImages(im, im2, best_transformation_matrix);
figure('Name','Best Stitched Image', 'FileName','BestStitchedImage.jpg');
imshow(bestStitchedImage);


% 1 (10 points) Hard Coding Point Correspondences
function [im, im2, points, points2] = point_correspondences(im, im2, points, points2)
    radius = 8;
    colors = {'red', 'yellow', 'blue', 'green'};
    for i = 1:size(points,1)
        im = insertShape(im, 'FilledCircle', [points(i,:) radius], 'Color', colors{i});
        im2 = insertShape(im2, 'FilledCircle', [points2(i,:) radius], 'Color', colors{i});
    end
end

%Author: Andrew Grier
%   - Finds the transformation matrix given a set of keypoints and
%     transformed points in the format:
%       [ x1 , y1 ]
%       [ x2 , y2 ]
%           ...
%       [ xn , yn ]
%     Outputs a transformation matrix specified in the format:
%       [ a , b , c ]
%       [ d , e , f ]
%       [ g , h , i ]
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

%Author: Andrew Grier
%   - Transforms multiple points coming in in the format:
%       [ x1 , y1 ]
%       [ x2 , y2 ]
%           ...
%       [ xn , yn ]
%     Using a transformation matrix specified in the format:
%       [ a , b , c ]
%       [ d , e , f ]
%       [ g , h , i ]
%     Outputs multiple points in the format:
%       [ x1' , y1' ]
%       [ x2' , y2' ]
%           ...
%       [ xn' , yn' ]
function transformedPoints = TransformMultiplePoints(points, transformationMatrix)
    transformedPoints = zeros(size(points, 1), size(points, 2));
    for pointNum = 1:size(points,1)
        transformedPoints(pointNum, :) = TransformPoint(points(pointNum,:), transformationMatrix);
    end
end

%Author: Andrew Grier
%   - Transforms a point coming in in the format:
%       [ x , y ]
%     Using a transformation matrix specified in the format:
%       [ a , b , c ]
%       [ d , e , f ]
%       [ g , h , i ]
%     Outputs a point in the format:
%       [ x , y ]
function transformedPoint = TransformPoint(point, transformationMatrix)
    %Assuming points are coming in as (x,y), homogenize the point.
    homogenousPoint = [
            point(1);
            point(2);
            1
        ];
    homogenousTransformedPoint = transformationMatrix * homogenousPoint;
    %Rounding here ensures the x and y are still in pixel amounts.
    transformedPoint = [round(homogenousTransformedPoint(1,:)/homogenousTransformedPoint(3,:)), round(homogenousTransformedPoint(2,:)/homogenousTransformedPoint(3,:))];
end

%Author: Andrew Grier
%   - Stitches two images together given a transformation matrix. This
%     function assumes that the transformed image is the one on the right
%     hand side. IMPORTANT!! THIS FUNCTION EXPECTS IMAGES IN THE UINT8
%     FORMAT FOR THEIR COLORS!!
function stitchedImage = StitchImages(baseImage, transformedImage, transformationMatrix)
    %Points are still in (x,y) format
    transformedImageCornerPoints = [
        1 1;                                                %Top Left
        1 size(transformedImage, 1);                        %Bottom Left
        size(transformedImage, 2) 1;                        %Top Right
        size(transformedImage, 2) size(transformedImage, 1) %Bottom Right
    ];

    %This part gets min and max values for x and y of the transformed image 
    %in the base image coordinate space to help determine what the size of 
    %the stitched image needs to be
    transformedImageTransformedCornerPoints = TransformMultiplePoints(transformedImageCornerPoints, inv(transformationMatrix));
    minTransformedY = min(transformedImageTransformedCornerPoints(1,2), transformedImageTransformedCornerPoints(3,2));
    maxTransformedY = max(transformedImageTransformedCornerPoints(2,2), transformedImageTransformedCornerPoints(4,2));
    minTransformedX = min(transformedImageTransformedCornerPoints(1,1), transformedImageTransformedCornerPoints(2,1));
    maxTransformedX = max(transformedImageTransformedCornerPoints(3,1), transformedImageTransformedCornerPoints(4,1));
    
    %This part gets min and max values for x and y of the base image 
    %to help determine what the size of the stitched image needs to be
    minBaseImageY = 0;
    maxBaseImageY = size(baseImage, 1);
    minBaseImageX = 0;
    maxBaseImageX = size(baseImage, 2);

    %This is used later for determining alpha to make edges in the images
    %appear less noticeable
    baseImageWidth = size(baseImage, 2);
    
    %Setup canvas for the stitched image
    minStitchedImageY = min(minBaseImageY, minTransformedY);
    maxStitchedImageY = max(maxBaseImageY, maxTransformedY);
    minStitchedImageX = min(minBaseImageX, minTransformedX);
    maxStitchedImageX = max(maxBaseImageX, maxTransformedX);
    stitchedImageHeight = maxStitchedImageY - minStitchedImageY;
    stitchedImageWidth = maxStitchedImageX - minStitchedImageX;

    %The output image will also be in uint8
    stitchedImage = uint8(zeros(stitchedImageHeight, stitchedImageWidth, 3));

    for y = 1:stitchedImageHeight
        %Simple Progress bar
        clc;
        disp("Stitching Progress = " + round((y/stitchedImageHeight)*100) + "%");
        for x = 1:stitchedImageWidth
            %If the transformed image made the canvas larger to the left of
            %the base image, we need to account for this when getting the
            %point from base image
            if (minTransformedX < 1)
                baseImageX = x + minTransformedX;
            else
                baseImageX = x;
            end
            
            %If the transformed image made the canvas taller than the top of
            %the base image, we need to account for this when getting the
            %point from base image
            if (minTransformedY < 1)
                baseImageY = y + minTransformedY;
            else
                baseImageY = y;
            end
            
            %Setup the two points to be in the same image space
            baseImagePoint = [baseImageX, baseImageY];
            transformedImagePoint = TransformPoint(baseImagePoint, transformationMatrix);
            transformedImageY = transformedImagePoint(2);
            transformedImageX = transformedImagePoint(1);
            
            %Check if the point we are looking for exists in both images,
            %if it doesn't we will need to set its color to black.
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
            %If the point exists in an overlap between the two images, we
            %need to blend their colors together, otherwise just use one
            %point or the other.
            if (pointExistsInBaseImage && pointExistsInTransformedImage)
                stitchedImage(y,x,:) = uint8(BlendPointColors(baseImagePointColor, transformedImagePointColor, (1 - baseImageX/baseImageWidth)));
            elseif (pointExistsInTransformedImage)
                stitchedImage(y,x,:) = uint8(transformedImagePointColor);
            else
                %Defaulting to baseImagePointColor also accounts for any
                %out of bounds points that are set to black above
                stitchedImage(y,x,:) = uint8(baseImagePointColor);
            end
        end
    end
end 

%Author: Andrew Grier
%   - Simple function that returns true if the point exists within the
%     bounds of the image. Needs both the point to be in the format:
%       [ x , y ]
function existsInImage = PointInImage (point, image)
    existsInImage = (point(1) > 0) && (point(1) <= size(image,2)) && (point(2) > 0) && (point(2) <= size(image,1));
end

%Author: Andrew Grier
%   - Given point colors in uint8, it will use the given alpha (between 0
%     and 1) to blend alpha% of the first color and 1-alpha% of the second
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

%Author: Andrew Grier
%   - Based on the given mode, find out if the point is the largest or
%     smallest out of all of its peers in a 3 by 3 by 3 area. The var mode
%     can be either "min" or "max" to determine if the point is a local
%     minimum or maximum respectively
%     Returns a boolean (0 or 1)
function isCurrentPointExtrema = IsCurrentPointExtrema(inputMatrix, mode)
    matrixHeight = size(inputMatrix, 1);
    matrixWidth = size(inputMatrix, 2);
    matrixDepth = size(inputMatrix, 3);
    currentValue = inputMatrix(2,2,2);
    isCurrentPointExtrema = false;
    for y = 1:matrixHeight
        for x = 1:matrixWidth
            for z = 1:matrixDepth
                % Skip the center point, since that is what we are 
                % comparing everything else to
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
    pruned_extrema_bits(1:border_distance, :) = 0;
    pruned_extrema_bits(end-border_distance+1:end, :) = 0;
    pruned_extrema_bits(:, 1:border_distance) = 0;
    pruned_extrema_bits(:, end-border_distance+1:end) = 0;

    std_dev_threshold = 0.03;
    patch_size = 9;
    patch_radius = floor(patch_size/2);

    [pruned_extrema_height, pruned_extrema_width] = find(pruned_extrema_bits == 1);
    for i = 1:length(pruned_extrema_height)
        patch = pruned_extrema_image(pruned_extrema_height(i) - patch_radius:pruned_extrema_height(i) + patch_radius,pruned_extrema_width(i) - patch_radius:pruned_extrema_width(i) + patch_radius);
        patch = double(patch(:));
        if std(patch(:)) < std_dev_threshold
            pruned_extrema_bits(pruned_extrema_height(i), pruned_extrema_width(i)) = 0;
        end
    end
    
    [pruned_extrema_height, pruned_extrema_width] = find(pruned_extrema_bits == 1);
    pruned_extrema_image = insertShape(baseImage, 'Circle', [pruned_extrema_width, pruned_extrema_height, ones(length(pruned_extrema_width), 1)*5], 'Color', 'red');
end

%Author: Andrew Grier
%   - Returns a specialized set, unionSet to be used with
%     keypoint_matching, this should be the union of the two input sets, 
%     minus any matches with distances that go above the threshold
function unionedSet = CreateUnionOfDescriptorSets(descriptorSet1, descriptorSet2, maxDistance)
    unionedSetTemp = [];
    unionedSetTempIndex = 1;
    %Get the union of the two sets
    for i = 1:size(descriptorSet1, 1)
        currentBestPairing = [];
        currentBestDistance = 5*maxDistance;
        for j = 1:size(descriptorSet2, 1)
            if i == descriptorSet2(j,2)
                %Prune descriptors with distances that are too far
                if (descriptorSet1(i,3) > maxDistance)
                    continue;
                end
                if (descriptorSet1(i,3) < currentBestDistance)
                    currentBestPairing = descriptorSet1(i, :);
                    currentBestDistance = descriptorSet1(i,3);
                end
            end
        end
        if (currentBestDistance <= maxDistance)
            unionedSetTemp(unionedSetTempIndex, :) = currentBestPairing;
            unionedSetTempIndex = unionedSetTempIndex + 1;
        end
    end
    
    unionedSet = [];
    unionedSetIndex = 1;
    %Remove duplicate entries from column 2
    for k = 1:size(unionedSetTemp, 1)
        currentPairedIndex = unionedSetTemp(k, 2);
        currentPairedDistance = unionedSetTemp(k,3);
        currentPairing = unionedSetTemp(k,:);
        foundBetterDuplicate = false;
        for m = 1:size(unionedSetTemp, 1)
            if (m == k)
                continue;
            end
            comparisonPairedIndex = unionedSetTemp(m,2);
            if (comparisonPairedIndex == currentPairedIndex)
                comparisonPairedDistance = unionedSetTemp(m,3);
                    if (comparisonPairedDistance < currentPairedDistance)
                        foundBetterDuplicate = true;
                    end
            end
        end
        if (foundBetterDuplicate)
            continue;
        else
            unionedSet(unionedSetIndex, :) = currentPairing;
            unionedSetIndex = unionedSetIndex + 1;
        end
    end
end

% 5 (10 points) Keypoint Description and Matching
function [canvas, keypoints1, keypoints2, C_union] = keypoint_matching(im, im2, remaining_extrema1, remaining_extrema2)
    %IMPORTANT: The indices for descriptors are the same as for their
    %keypoints locations. So we can use these later when drawing lines
    [keypoints1, descriptors1] = extract_descriptors(im, remaining_extrema1);
    [keypoints2, descriptors2] = extract_descriptors(im2, remaining_extrema2);

    C_1 = match(descriptors1, descriptors2);
    C_2 = match(descriptors2, descriptors1);

    %Indices from keypoints1 are on the left
    %Indices from keypoints2 are in the middle
    %Distances are on the right
    %   Feel free to change the third input to a threshold of your choice!
    C_union = CreateUnionOfDescriptorSets(C_1, C_2, 300);
    
    %Remember to look up how drawing the lines will change based on the
    %wider canvas with both images...
    canvas = draw_matches(im, im2, keypoints1, keypoints2, C_union);
end

function [keypoints, descriptors] = extract_descriptors(im, remaining_extrema)
    [extrema_height, extrema_width] = find(remaining_extrema == 1);
    keypoints = [extrema_height, extrema_width];
    extrema_length = length(extrema_height);
    descriptors = zeros(extrema_length, 243);
    for i = 1:extrema_length
        patch = im(extrema_height(i)-4:extrema_height(i)+4, extrema_width(i)-4:extrema_width(i)+4, :);
        flatten = reshape(patch, [1, 243]);
        descriptors(i,:) = flatten;
    end
end

%Author: Andrew Grier
%   - Takes in the two sets of descriptors and matches the closest ones
%     using the indices of the two lists.
%     Returns in the format:
%       [ descriptor1Index1 , descriptor2IndexOfBestMatch1 , distanceForThisMatch1 ]
%       [ descriptor1Index2 , descriptor2IndexOfBestMatch2 , distanceForThisMatch2 ]
%                                           ...
%       [ descriptor1IndexN , descriptor2IndexOfBestMatchN , distanceForThisMatchN ]
function matches = match(descriptors1, descriptors2)
    descriptors_height = size(descriptors1, 1);
    descriptors2_height = size(descriptors2, 1);
    % The index of the source point is on the left, the index of the
    % matched point is in the middle,
    % distance is on the right.

    matches = zeros(descriptors_height, 3);
    
    for i = 1:descriptors_height
        keypointDescriptorSource = descriptors1(i,:);
        % If anything still has an index of -1 by the end, something went
        % wrong!
        currentBestMatchIndex = -1;
        currentBestMatchDistance = 999999;
        for j = 1:descriptors2_height
            keypointDescriptorComparison = descriptors2(j,:);
            distance = sqrt((keypointDescriptorSource - keypointDescriptorComparison) * (keypointDescriptorSource - keypointDescriptorComparison).');
            if (distance <= currentBestMatchDistance)
                currentBestMatchIndex = j;
                currentBestMatchDistance = distance;
            end
        end
        matches(i,1) = i;
        matches(i,2) = currentBestMatchIndex;
        matches(i,3) = currentBestMatchDistance;
    end
end

function canvas = draw_matches(im, im2, keypoints1, keypoints2, C_union)
    canvas_height = max(size(im, 1), size(im2, 1));
    canvas_width = size(im, 2) + size(im2, 2);
    canvas = zeros(canvas_height, canvas_width, size(im, 3), 'like', im);
    canvas(1:size(im, 1), 1:size(im, 2), :) = im;
    canvas(1:size(im2, 1), size(im, 2)+1:end, :) = im2;

    keypoints2(:,2) = keypoints2(:,2) + size(im, 2);
    for i = 1:size(C_union, 1)
        start_point = keypoints1(C_union(i, 1), :);
        end_point = keypoints2(C_union(i, 2), :);
        canvas = insertShape(canvas, 'Line', [[start_point(2), start_point(1)], [end_point(2), end_point(1)]], 'Color', 'red', 'LineWidth', 2);
    end
end

% 6 (10 points) Find the Transformation Matrix via RANSAC and Stitch
function [canvas, best_transformation_matrix] = auto_stitch(im, im2, keypoints1, keypoints2, C_union)
    best_transformation_matrix = [];
    threshold_distance = 20;
    best_close_distances = 0;
    im_size = size(im);
    im_size2 = size(im2);
    best_random_extrema1 = zeros(size(im_size));
    best_random_extrema2 = zeros(size(im_size2));

    experiments = 100;
    number_of_correspondences = 4;
    for i = 1:experiments
        random_indices = randperm(size(C_union, 1), number_of_correspondences);
        random_keypoints1 = keypoints1(C_union(random_indices,1), :);
        random_keypoints2 = keypoints2(C_union(random_indices,2), :);
        random_extrema1 = zeros(im_size);
        random_extrema2 = zeros(im_size2);

        current_transformation_matrix = FindTransformationMatrixWithPoints(flip(random_keypoints1,1), flip(random_keypoints2,1));

        current_close_distances = 0;
        for j = 1:size(C_union,1)
            current_keypoint1 = keypoints1(C_union(j,1),:);
            current_keypoint2 = keypoints2(C_union(j,2),:);
            transformed_point = TransformPoint(current_keypoint1, current_transformation_matrix);
            distance = sqrt(((transformed_point(1) - current_keypoint2(1))^2) + ((transformed_point(2) - current_keypoint2(2))^2));
            if distance <= threshold_distance
                current_close_distances = current_close_distances + 1;
            end
        end

        if isempty(best_transformation_matrix) || current_close_distances > best_close_distances
            best_transformation_matrix = current_transformation_matrix;
            for j = 1:number_of_correspondences
                random_extrema1(random_keypoints1(j,1), random_keypoints1(j,2)) = 1;
                random_extrema2(random_keypoints2(j,1), random_keypoints2(j,2)) = 1;
            end
            best_random_extrema1 = random_extrema1;
            best_random_extrema2 = random_extrema2;
        end
    end
    canvas = keypoint_matching(im, im2, best_random_extrema1, best_random_extrema2);
end
