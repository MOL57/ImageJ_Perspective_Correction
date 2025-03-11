# ImageJ_Perspective_Correction
**ImageJ macro to correct/rectify the perspective of an image.**

Requires ImageJ version 1.51r  or higher.

Corrects the perspective by solving a homographic transformation
given the coordinates of 4 points on the input image
and the coordinates of the 4 corresponding points on the output image.

![Correcting perspective process](images/fig1.jpg)

The homography is solved using the method described in:\
https://math.stackexchange.com/questions/494238/how-to-compute-homography-matrix-h-from-corresponding-points-2d-2d-planar-homog

The system of linear equations is solved using the Gauss-Jordan method as described in Numerical recipes in C, see:\
http://s3.amazonaws.com/nrbook.com/book_C210.html \
and:\
https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/InClassExamples/NR3/legacy/nr2/C_211/recipes/gaussj.c

