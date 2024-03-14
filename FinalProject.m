% PLEASE CHANGE THIS CODE IT IS UGLY

clear variables;
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

%%%%%%%%%%%%%%
% QUESTION 1 %
%%%%%%%%%%%%%%
figure('Name','Point Correspondence Side-by-Side', 'FileName','PointCorrespondence.jpg');
imshowpair(im,im2,'montage');

%%%%%%%%%%%%%%
% QUESTION 2 %
%%%%%%%%%%%%%%
transformationMatrix = FindTransformationMatrixWithPoints(points, points2);


% 1 (10 points) Hard Coding Point Correspondences
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

% Get the transformation matrix using points
function homogenousMatrix = FindTransformationMatrixWithPoints(points, points2)
    %Note: Points are still in (x,y) notation
    
    homogenousMatrix = zeros(3,3);
    A = [];
    B = [];
    rowNum = 1;
    for pointNum = 1:4
        A(rowNum, :) = [points(pointNum, 1) points(pointNum, 2) 1 0 0 0 0 0 0];
        B(rowNum, :) = [0 0 0 0 0 0 (points(pointNum, 1) * points2(pointNum, 1)) (points(pointNum, 2) * points2(pointNum, 1)) points2(pointNum, 1)];
        rowNum = rowNum + 1;
        A(rowNum, :) = [0 0 0 points(pointNum, 1) points(pointNum, 2) 1 0 0 0];
        B(rowNum, :) = [0 0 0 0 0 0 (points(pointNum, 1) * points2(pointNum, 2)) (points(pointNum, 2) * points2(pointNum, 2)) points2(pointNum, 2)];
        rowNum = rowNum + 1;
    end
    At = A.';
    Bt = B.';
    [U, S, V] = svd((At * A) + (Bt * B) - (Bt * A) - (At * B));
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
    testPoint = points(testPointNum,:);
    testPointPrime = points2(testPointNum, :);
    transformTestPoint = TransformPoint(testPoint, homogenousMatrix);
    isTransformationCorrect = isequal(testPointPrime, transformTestPoint);
    if (~isTransformationCorrect) 
        disp("Something went wrong finding the transformation matrix!");
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