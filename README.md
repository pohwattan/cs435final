Andrew Grier
Nick Pohwat

Final Project
3/21/24

1. Features of your program:

	All parts of the code are commented with their respective questions.

	a. Hard Coding Point Correspondences
		- Will take in two images and two arrays of four points and add marks.
		- Returns: two images marked with point correspondences.
		- Displays the two images marked up.

	b. Compute Transformation Matrix, Project, and Blend!
		- Takes two arrays of four points.
		- Returns: A transformation matrix for the points.
		- Displays the stitched images

	c. Create Scale-Space Image Pyramids
		- Takes in two images.
		- Returns: Two image pyramids and two Difference of Gaussian pyramids.
		- Displays the two image pyramids.

	d. Finding the Local Maximas
		- Takes in a Difference of Gaussian pyramid and an original grayscale image.
		- Returns: An all extrema image, a pruned extrema image, and the remaining extrema.
		- Displays two figures, each with the all extrema image and the pruned extrema image.

	e. Keypoint Description and Matching
		- Takes the two images and the remaining extrema for each image.
		- Returns: The keypoint matched image, the keypoints, and the C_union.
		- Displays the keypoint matched image.

	f. Find the Transformation Matrix via RANSAC and Stitch
		- Takes the images, the keypoints, and the C_union.
		- Returns: the best keypoints matched images and the best transformation matrix.
		- Displays the best keypoints matched images and the stitched image.

2. Name of your entry-point script:

	FinalProject.m

3. Any useful instructions to run your script:

	a. Have any images within the same directory level as FinalProject.m to run as is with our images.
	b. Image names should be "first.jpg" and "second.jpg"
	c. Images are read in on lines 4 and 5.