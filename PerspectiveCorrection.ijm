////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// CORRECTING THE PERSPECTIVE OF AN IMAGE
// by means of solving a homographic transformation
// based on 4 reference points on the input image
// which correspond, on the (perspective-corrected) output image,
// to the the vertices of a square or rectangle.
//
// Optionally the 4 reference points on the output image
// can be placed on general positions (instead of on the vertices of a rectangle/square)
// thus allowing for an image warping instead of perspective correction
//
// PerspectiveCorrection.ijm
// 2025.02.25
// (c) Xavier Fernandez
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// MAIN BODY OF THE MACRO
/////////////////////////////////////////////////////////////////////////////////

var VERBOSE = 0;		// print progress on the log window

requires("1.51r");

if(VERBOSE)	print("CORRECTING THE PERSPECTIVE OF AN IMAGE\n");

// INPUT IMAGE AND METHOD TO SELECT THE REFERENCE POINTS

openImages = getList("image.titles");
if(openImages.length==0)	openImages[0] = "-- no images --";

choices_input = newArray("graphically", "numerically");
choices_output = newArray("graphically (points on vertices of a rectangle/square)", "graphically (points on arbitrary positions)", "numerically");

Dialog.createNonBlocking("Perspective correction");
Dialog.addMessage("CORRECTING THE PERSPECTIVE OF AN IMAGE\n \n" +
			"The input image (to be corrected) must contain 4 reference points\n" +
			"which correspond to 4 reference points on the output image (corrected)");
Dialog.addChoice("Image to be corrected:", openImages);	
Dialog.addRadioButtonGroup("Method to enter reference points on INPUT image:", choices_input, 2, 1, "graphically");
Dialog.addRadioButtonGroup("Method to enter reference points on OUTPUT image:", choices_output, 3, 1, "graphically (points on vertices of a rectangle/square)");
Dialog.show;
imgName = Dialog.getChoice();
inputPointsMethod = Dialog.getRadioButton;
outputPointsMethod = Dialog.getRadioButton;

if( ! isOpen(imgName) )	exit("ERROR: Image not found");
selectImage(imgName);

if(VERBOSE)
{
	print("Input image:",imgName);
	print("Selection method for points on input image:", inputPointsMethod);
	print("Selection method for points on output image:", outputPointsMethod);
	print("");
}


// SELECTING THE REFERENCE POINTS ON THE INPUT IMAGE

if(inputPointsMethod == "graphically")	// selecting points graphically
{
	setTool("polygon");

	Dialog.createNonBlocking("Points on input image");
	Dialog.addMessage("Draw on the image a polygon\n" +
				"whose 4 vertices mark the positions of the references\n" +
				"starting from the top-left reference point and going clockwise:");
	Dialog.addMessage("start > > > > > > >1\n"+
				 "  ^                            v\n"+
				 "  ^                            v\n"+
				 "  ^                            v\n"+
				 "  3 < < < < < < < < 2"); 
	Dialog.addMessage("Press OK when done\n ");
	Dialog.show;

	selectImage(imgName);
	type = selectionType;
	//print("Selection type is:", type);
	if( type<0 )	exit("ERROR: Image contains no selection");
	if( type!=2 )	exit("ERROR: Selection must be a polygon");

	getSelectionCoordinates(U, V);
	n = U.length;
	if( n!=4 )	exit("ERROR: Selection must have 4 points");
}
else	// entering points numerically
{
	U = newArray(4);
	V = newArray(4);
	n = 4;

	Dialog.createNonBlocking("Points on input image");
	Dialog.addMessage("Enter the coordinates of the reference points\non the INPUT image");
	Dialog.addMessage("(U0,V0) ---------------- (U1,V1)\n"+
				"      |                                    |      \n" +
				"      |                                    |      \n" +
				"      |                                    |      \n" +
				"(U3,V3) ---------------- (U2,V2)\n");		
	Dialog.addNumber("U0", 0); Dialog.addToSameRow(); Dialog.addNumber("V0", 0);
	Dialog.addNumber("U1", 0); Dialog.addToSameRow(); Dialog.addNumber("V1", 0);
	Dialog.addNumber("U2", 0); Dialog.addToSameRow(); Dialog.addNumber("V2", 0);
	Dialog.addNumber("U3", 0); Dialog.addToSameRow(); Dialog.addNumber("V3", 0);
	Dialog.show;
	
	for(i=0; i<4; i++)
	{
		U[i] = Dialog.getNumber();
		V[i] = Dialog.getNumber();
	}
}
	
if(VERBOSE)
{
	print("Coordinates of the reference points\non the input (uncorrected) image:");
	print("#      U     V");
	for(i=0;i<n;i++)	print(i, "  ", U[i], "  ", V[i]);
}

// SELECTING THE REFERENCE POINTS ON THE OUTPUT IMAGE

if(outputPointsMethod == "numerically")	// entering points numerically
{
	X = newArray(4);
	Y = newArray(4);
	n = 4;

	Dialog.createNonBlocking("Points on output image");
	Dialog.addMessage("Enter the coordinates of the reference points\non the OUTPUT image");
	Dialog.addMessage("(X0,Y0) ---------------- (X1,Y1)\n"+
				"      |                                    |      \n" +
				"      |                                    |      \n" +
				"      |                                    |      \n" +
				"(X3,Y3) ---------------- (X2,Y2)\n");	
	Dialog.addNumber("X0", 0); Dialog.addToSameRow(); Dialog.addNumber("Y0", 0);
	Dialog.addNumber("X1", 0); Dialog.addToSameRow(); Dialog.addNumber("Y1", 0);
	Dialog.addNumber("X2", 0); Dialog.addToSameRow(); Dialog.addNumber("Y2", 0);
	Dialog.addNumber("X3", 0); Dialog.addToSameRow(); Dialog.addNumber("Y3", 0);
	Dialog.show;
	
	for(i=0; i<4; i++)
	{
		X[i] = Dialog.getNumber();
		Y[i] = Dialog.getNumber();
	}
}
else if(outputPointsMethod == "graphically (points on vertices of a rectangle/square)")	// selecting points graphically on a rectangle/square
{
	selectImage(imgName);
	run("Select None");

	setTool("rectangle");

	Dialog.createNonBlocking("Points on output image");
	Dialog.addMessage("Draw on the image a square/rectangle\n" +
				"whose 4 vertices mark the goal positions of the references\n" +
				"once the image perspective had been corrected");
	Dialog.addMessage("Press and hold <Shift> while drawing to force a square\nPress OK when done\n ");				
	Dialog.addMessage("start > > > > > > >1\n"+
				 "  ^                            v\n"+
				 "  ^                            v\n"+
				 "  ^                            v\n"+
				 "  3 < < < < < < < < 2"); 
	Dialog.show;

	type = selectionType;
	//print("Selection type is:", type);
	if( type<0 )	exit("ERROR: Image contains no selection");
	if( type!=0 )	exit("ERROR: Selection must be a rectangle/square");

	getSelectionCoordinates(X, Y);
	n = X.length;
	if( n!=4 )	exit("ERROR: Selection must have 4 points");
}
else 	// selecting points graphically on a polygon
{
	selectImage(imgName);
	run("Select None");

	setTool("polygon");

	Dialog.createNonBlocking("Points on output image");
	Dialog.addMessage("Draw on the image a polygon\n" +
				"whose 4 vertices mark the goal positions of the references\n" +
				"(once the image perspective had been corrected)\n" +
				"starting from the top-left reference point and going clockwise:");				
	Dialog.addMessage("start > > > > > > >1\n"+
				 "  ^                            v\n"+
				 "  ^                            v\n"+
				 "  ^                            v\n"+
				 "  3 < < < < < < < < 2"); 
	Dialog.show;

	type = selectionType;
	//print("Selection type is:", type);
	if( type<0 )	exit("ERROR: Image contains no selection");
	if( type!=2 )	exit("ERROR: Selection must be a polygon");

	getSelectionCoordinates(X, Y);
	n = X.length;
	if( n!=4 )	exit("ERROR: Selection must have 4 points");
}

if(VERBOSE)
{
	print("\nCoordinates of the reference points\non the corrected image:");
	print("#      X     Y");
	for(i=0;i<n;i++)	print(i, "  ", X[i], "  ", Y[i]);
}

// Message box while processing

newImage("Processing", "RGB", 270, 70, 1);
setColor(240,240,240);
fill();
setColor(0,0,0);
setFont("SansSerif", 15);
drawString("Wait while correcting perspective ...", 20, 45);
setLocation(screenWidth/3, screenHeight/3);

// SOLVING THE HOMOGRAPHY MATRIX

H = newArray( 1, 0, 0,		// 3X3 matrix holding the coefficients of the transformation of the input image to the output image
		 0, 1, 0,		// matrix expressed as a unidimensional array in row order: h00, h01, h02, h10, h11, h12, h20, h21, h22
		 0, 0 ,1); 

solveHomography4points(X, Y, U, V, H);

if(VERBOSE)
{
	print("Matrix \"H\" of coefficients of the homographic transformation:");
	print (H[0*3+0], ", ", H[0*3+1], ", ", H[0*3+2]);
	print (H[1*3+0], ", ", H[1*3+1], ", ", H[1*3+2]);
	print (H[2*3+0], ", ", H[2*3+1], ", ", H[2*3+2]);
	print("");
}

// ERROR ESTIMATION:
// applying the transformation matrix to the X, Y coordinates
// should yield the U, V coordinates

U_ = newArray( 0, 0, 0, 0);	// computed output coordinates
V_ = newArray(0, 0, 0, 0);	

for( point=0; point<4; point++)
{
	U_[point] = H[ 0*3 + 0] * X[point] + H[ 0*3 + 1] * Y[point] + H[ 0*3 + 2];
	V_[point] = H[ 1*3 + 0] * X[point] + H[ 1*3 + 1] * Y[point] + H[ 1*3 + 2];
		 Z = H[ 2*3 + 0] * X[point] + H[ 2*3 + 1] * Y[point] + H[ 2*3 + 2];

	U_[point] = U_[point] / Z;
	V_[point] = V_[point] / Z;
}

if(VERBOSE)
{
	print("Error estimation:");

	for(point=0; point<4; point++)
	{
		print( "point", point, ":", "X:", X[point], "=> U_:", U_[point], " ( U:", U[point], ", error:", U_[point] - U[point],")" );
		print( "point", point, ":", "Y:", Y[point], "=> V_:", V_[point], " ( V:", V[point], ", error:", V_[point] - V[point],")" );
	}
	print("");
}

// APPLYING THE HOMOGRAPHY TO THE INPUT IMAGE

// Building the macro expression which maps the coordinates on the input image
// onto the coordinates on the output image

text = macroExprHomography(H);

if(VERBOSE)
{
	print("Macro expression:");
	print(text);
	print("");
}

// Apply transform to the image

setBatchMode(true);

selectImage(imgName);
run("Select None");
correctedName = "CORRECTED_" + imgName;
run("Duplicate...", "title=[" + correctedName + "]");

selectImage(correctedName);
run("Macro...", text);

setBatchMode(false);

// SHOWING RESULTS

selectImage(imgName);
makeSelection("polygon", U, V);

selectImage(correctedName);
makeSelection("polygon", X, Y);

if(VERBOSE)	print("IMAGE PERSPECTIVE SUCCESSFULLY CORRECTED");

selectImage("Processing");
close();
waitForUser("Perspective successfully corrected!");

close("Log");
exit;

//////////////////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION TO SOLVE THE HOMOGRAPHY MATRIX
/////////////////////////////////////////////////////////////////////////////////////////////////

// Solves the homography matrix of a perspective projection
// given the coordinates of 4 points on the first projection
// and the coordinates of the 4 corresponding points on the second projection

//
// SOURCE:
// https://towardsdatascience.com/estimating-a-homography-matrix-522c70ec4b2c
// https://www.cs.umd.edu/class/fall2019/cmsc426-0201/files/16_Homography-estimation-and-decomposition.pdf
//
// Arguments:
//	X:	Array holding the X coordinates of the 4 points on the first projection
//	Y:	Array holding the Y coordinates of the 4 points on the first projection
//	U:	Array holding the X coordinates of the 4 points on the second projection
//	V:	Array holding the Y coordinates of the 4 points on the second projection
//	H:	3x3 matrix expressed as a a unidimensional array of 9 elements in row order,
//		that will contain, on return, the coefficients of the homography transformation 

function solveHomography4points(X, Y, U, V, H)
{

	a = newArray( 8*8 );		// 8x8 matrix containing the coefficients of the system of equations to be solved
					// expressed as a unidimensional array in row order
			 
	b = newArray( 8 );		// vector containing the 8 independent coefficients of the system of equations to be solved

	// Filling values of "a" matrix

	for( point = 0; point<4; point++)
	{
		a[ point*8 + 0 ] = X[ point ];
		a[ point*8 + 1 ] = Y[ point ];
		a[ point*8 + 2 ] = 1;
		a[ point*8 + 3 ] = 0;
		a[ point*8 + 4 ] = 0;
		a[ point*8 + 5 ] = 0;
		a[ point*8 + 6 ] = - X[ point ] * U[ point ];
		a[ point*8 + 7 ] = - Y[ point ] * U[ point ];

		a[ (4+point)*8 + 0 ] =0;
		a[ (4+point)*8 + 1 ] =0;
		a[ (4+point)*8 + 2 ] =0;
		a[ (4+point)*8 + 3 ] =X[ point ];
		a[ (4+point)*8 + 4 ] =Y[ point ];
		a[ (4+point)*8 + 5 ] =1;
		a[ (4+point)*8 + 6 ] = - X[ point ] * V[ point ];
		a[ (4+point)*8 + 7 ] = - Y[ point ] * V[ point ];
	}

	// Filling values of "b" vector

	for( point = 0; point<4; point++)
	{
		b[ point ] = U[ point ];
		b[ point + 4 ] = V[ point ];
	}

	if(VERBOSE)
	{
		print("");	
		print("Matrix \"a\" of coefficients of the system of equations:");
		for( row=0; row<8; row++)
			print(a[row*8+0],", ", a[row*8+1],", ", a[row*8+2],", ", a[row*8+3],", ", a[row*8+4],", ", a[row*8+5],", ", a[row*8+6],", ", a[row*8+7]);
		
		print("");
		print("Vector \"b\" of independent values of the system of equations:");
		print(b[0],", ", b[1],", ", b[2],", ", b[3],", ", b[4],", ", b[5],", ", b[6],", ", b[7]);
		print("");
	}

	// Solving the system of equations

	solveSystLinEq(a, 8, b);

	// Building H matrix from the returned values in the "b" vector

	for( i=0; i<8; i++)	H[i] = b[i];
	H[8] = 1;

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION TO SOLVE A SYSTEM OF LINEAL EQUATIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Solving a system of linear equations by means of Gauss-Jordan method with pivot
//
// SOURCE:
// Numerical recipes in C
// http://s3.amazonaws.com/nrbook.com/book_C210.html
// https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/InClassExamples/NR3/legacy/nr2/C_211/recipes/gaussj.c
//
//	a00 * x0 + a01 * x1 + a02 * x2 + ... + a0(n-1) * x(n-1) = b0
//	a10 * x0 + a11 * x1 + a12 * x2 + ... + a1(n-1) * x(n-1) = b1
//			...					...
//	a(n-1)0 * x0 + a(n-1)1 * x1 + ...+a(n-1)(n-1) * x(n-1) = b(n-1)
//
// Arguments:
//	a:		Matrix containing the coefficients aij of the equations
//			expressed as a unidimensional array in row order:
//			a00, a01,..a0(n-1), a10, a11, ...
//	n:		Number of equations (and unknowns)
//	b:		Array containing the independent coefficients: b0, b1... b(n-1)
//			On return it will contain the solutions: x0, x1, ... x(n-1)


function solveSystLinEq(a, n, b)
{
	N_MAX = 10;		// Maximum number of equations and unknowns

	if( n > N_MAX)	exit("ERROR: too much equations");

	i = 0; icol = 0; irow = 0; j = 0; k = 0; l = 0; ll = 0;	
	big = 0.0; dum = 0.0; pivinv = 0.0; temp = 0.0;	
	
	indxc = newArray(N_MAX);	// The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting
	indxr = newArray(N_MAX);
	ipiv = newArray(N_MAX);

	for( j = 0; j < n; j++)	ipiv[ j ] = 0;	

	for( i = 0 ;i < n; i++)			//This is the main loop over the columns to be reduced 
	{			
		big = 0.0;
		for( j=0; j < n; j++)		// This is the outer loop of the search for a pivot element 
			if( ipiv[j] != 1)
				for( k=0; k < n; k++)
				{
					if( ipiv[k] == 0)
					{
						if( abs( a[ j*n + k ] ) >= big)
						{
							big = abs( a[ j*n + k ]);
							irow = j;
							icol = k;
						}
					}
				}
		ipiv[icol] = ipiv[icol] +1;

				
		// We now have the pivot element, so we interchange rows, if needed, to put the pivot
		//element on the diagonal. The columns are not physically interchanged, only relabeled:
		//indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
		//indxr[i] is the row in which that pivot element was originally located. If indxr[i] _=
		//indxc[i] there is an implied column interchange. With this form of bookkeeping, the
		//solution b’s will end up in the correct order, and the inverse matrix will be scrambled
		//by columns. 

		if( irow != icol )
		{
			for( l=0; l<n; l++)
			{
				temp = a[ irow*n + l ];
				a[ irow*n + l ] = a [icol*n + l ];
				a [icol*n + l ] = temp;
			}

			temp = b[ irow ];
			b[ irow ] = b[ icol ];
			b[ icol ] = temp;
	
		}

		indxr[i]=irow;							//We are now ready to divide the pivot row by the
		indxc[i]=icol;							//pivot element, located at irow and icol.

		//print(irow, icol);

		if( a[ icol*n + icol ] == 0.0)	exit("ERROR: Singular Matrix");
		
		pivinv = 1.0 / a[ icol*n + icol ];
		a[ icol*n + icol ]=1.0;
		for( l = 0; l < n; l++) 	a[ icol*n + l ] = a[ icol*n + l ] * pivinv;
		b[ icol ] = b[ icol ] * pivinv;

		for( ll = 0; ll < n; ll++ )
			if( ll != icol )
			{
				dum = a[ ll*n + icol ];
				a[ ll*n + icol ]=0.0;
				for( l = 0; l < n; l++)	a[ ll*n + l ] = a[ ll*n + l ] - a[ icol*n + l ] * dum;
				b[ ll ] = b[ ll ] - b[ icol ] * dum;
			}
	}
	
	//This is the end of the main loop over columns of the reduction. It only remains to unscramble
	//the solution in view of the column interchanges. We do this  by interchanging pairs of
	//columns in the reverse order that the permutation was built up

	for( l = n-1; l >= 0; l--)
	{
		if( indxr[ l ] != indxc[ l ] )
			for( k = 0; k < n; k++ )
			{
				temp = a[ k*n + indxr[ l ] ];
				a[ k*n + indxr[ l ] ] = a[ k*n + indxc[ l ] ];
				a[ k*n + indxc[ l ] ] = temp;
			}
	}

	/* And we are done. */
}

/////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION TO BUILD THE MACRO EXPRESSION 
// TO APPLY HOMOGRAPHY TO AN IMAGE
/////////////////////////////////////////////////////////////////////////////////////////////
//
// Builds the text expression of the ImageJ macro to apply homography to an image
// Example:
//
//	run("Macro...", "code=[EXPRESSION]"); 
//
// The expression maps the value of each pixel in the original image
// to the same pixel value of the corresponding pixel in the transformed image
//
// Arguments:
//	H:		3x3 matrix expressed as a unidimensional array in row order
//			containing the coefficients of the homography transform
//
// Returns:		String with the text of the macro expression
//

function macroExprHomography(H)
{
	u_ = "( (" + H[ 0*3 + 0] +") * x + (" + H[ 0*3 + 1] +") * y + (" + H[ 0*3 + 2] +") ) / ( (" + H[ 2*3 + 0] +") * x + (" + H[ 2*3 + 1] +") * y + (" + H[ 2*3 + 2] +") )";
	v_ = "( (" + H[ 1*3 + 0] +") * x + (" + H[ 1*3 + 1] +") * y + (" + H[ 1*3 + 2] +") ) / ( (" + H[ 2*3 + 0] +") * x + (" + H[ 2*3 + 1] +") * y + (" + H[ 2*3 + 2] +") )";

	expr = "code=[v=getPixel(" + u_ + ", " + v_ + ")]";

	return expr;
}





