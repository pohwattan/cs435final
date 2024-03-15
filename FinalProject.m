% PLEASE CHANGE THIS CODE IT IS UGLY

clear variables;
close all;

im = imread("first.jpg");
im2 = imread("second.jpg");

% PREPROCESSING - Rotate images to vertical orientation
im = imrotate(im, -90);
im2 = imrotate(im2, -90);

[imMarked, im2Marked, points, points2] = point_correspondences(im, im2);

% figure(1);
% imshow(im);
% figure(2);
% imshow(im2);

%%%%%%%%%%%%%%
% QUESTION 1 %
%%%%%%%%%%%%%%
figure('Name','Point Correspondence Side-by-Side', 'FileName','PointCorrespondence.jpg');
imshowpair(imMarked,im2Marked,'montage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 2 by Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transformationMatrix = FindTransformationMatrixWithPoints(points, points2);
stitchedImage = StitchImages(im, im2, transformationMatrix);
figure('Name','Stitched Image', 'FileName','StitchedImage.jpg');
imshow(stitchedImage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 6 by Andrew Grier %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1 (10 points) Hard Coding Point Correspondences
function [im, im2, points, points2] = point_correspondences(im, im2)
    points = [2188 943; 1193 2036; 1125 3592; 2660 3585];
    points2 = [1541 1024; 463 2026; 257 3668; 1774 3577];
    radius = 20;
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