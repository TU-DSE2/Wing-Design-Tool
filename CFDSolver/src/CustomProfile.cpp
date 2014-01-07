#include "../include/CustomProfile.h"

/*CustomProfile::CustomProfile(mat Gamma_I_in, mat n_grid_in, mat interior_in, mat ex_loc_in, mat ex_n_in, vec ex_ds_in) {
    Gamma_I = Gamma_I_in;
    n_grid = n_grid_in;
    interior = interior_in;
    ex_loc = ex_loc_in;
    ex_n = ex_n_in;
    ex_ds = ex_ds_in;
}*/

CustomProfile::CustomProfile() {
    //Variables
    float gridsize = 32;
    float gridsizex=gridsize, gridsizey=gridsize;
    float alpha = 15, scale = (gridsize/5)*3;

    // Opening file using the function open_file, returning int Nrows
    // and matrix A, containing the coördinates
    string myfilename="374.dat";
    ifstream ifile(myfilename.c_str());
    if(ifile) {
        pair <mat, int> A_Nrows;
        A_Nrows=open_file(myfilename);          // Opening the file
        mat A = A_Nrows.first;                  // Seperating pair
        int Nrows = A_Nrows.second;

        // Define grid
        mat Grid = zeros(gridsizex,gridsizey);
        // Scaling the coods from A down and setting mean of y to half the grid
        A = scalecoodsdown(A);
        // Transform coordinates with the chosen angle of attack
        rotatecoods(alpha, A, Nrows);                    // Returns A, with rotated coods, not bigger then one
        // Scale the airfoil to the grid.
        scalecoods(scale, gridsize, A);                  // Returns A, with scaled coods, hopefully not bigger then the grid


        // Find the exact boundary
        mat ExactBndry = ExactBoundary(A, Nrows);

        // Get the original normal unitvectors from the airfoil itself
        mat normals = createoriginalnormals(1,A,Nrows);

        // Create the grid boundary matrices, containing the information in matrix format
        cube GridBndryCube = CreateGridBndryCube(A, normals, Nrows, gridsizex, gridsizey);
        // Find amount of gridpoints
        mat matGamma = GridBndryCube.slice(0);
        int nGamma = accu(matGamma);

        // Creating the Grid Bndry from the cube. Ordering them as well.
        mat GridBndry = CreateGridBndry(GridBndryCube,nGamma,gridsizex,gridsizey);

        // Find interior boundary points Omega,using pnpoly from the web.
        interior = FindOmega(GridBndry, nGamma, gridsizex, gridsizey);

        // Redefine matrices for saving
        Gamma_I = GridBndry.cols(0,1);
        n_grid = GridBndry.cols(2,3);

        ex_loc = ExactBndry.cols(0, 1);
        ex_ds = ExactBndry.col(2);
        ex_n = ExactBndry.cols(3, 4);
    } else {
        cout << "Error: file not found!";
    }
}

pair<mat, int> CustomProfile::open_file(string myfilename)
{
    // Define initial variables and read in file
        string line, line1, line2;
        int i = 0, Nrows = 0, number_of_lines =0;
        ifstream myfile (myfilename.c_str());

    // Find number of lines and define the A Matrix
        while (getline(myfile, line))
            {
                if ( line.empty() )
                {
                    continue;
                }
                number_of_lines++;

        }
        myfile.clear();                 // clears the eof bits and so on.
        Nrows=number_of_lines-1;
        mat A = zeros(Nrows,2);

    // Find size of the file and reset eof flags and set position
    // in the stream
        myfile.seekg (0, myfile.end);
        //int length = myfile.tellg(); (UNUSED)
        myfile.seekg (0, myfile.beg);

    // If myfile is open, it reads out the data
        if(myfile.is_open() && myfile.good())
            {
//                while(! myfile.eof() )
                while( i < number_of_lines)
                    {
                        if (i == 0)         //This loop reads the first line
                            {
                                getline(myfile, line);
                            }
                        else                //Here the data is read and put into matrix
                            {
                                getline(myfile, line1, ',');        // Reads out x
                                stringstream(line1) >> A(i-1,0);
                                getline(myfile, line2, '\n');
                                stringstream(line2) >> A(i-1,1) ;   // Reads out y
                            }
                        i++;
                    }
                myfile.close();         // Close the file
            }
    return make_pair(A,Nrows);
}

////////////////////////////////////////////////////////////////////////////////
// Rotating the input coordinates with a angle alpha clockwise.
// It is assumed that the input coods for x are from 0 to 1 and
// for y range from -0.5 to 0.5 (or that its mean is through the x axis)
void CustomProfile::rotatecoods(float alpha, mat& A, int Nrows)
{

    mat Atmp = zeros(Nrows,1);
    alpha =alpha/360*2*M_PI;
    if (alpha != 0)
    {
        Atmp.col(0)= ( cos(alpha) * (A.col(0)-0.5) ) + ( sin(alpha) * A.col(1) );
        A.col(1)=    ( -sin(alpha)* (A.col(0)-0.5) ) + ( cos(alpha) * A.col(1) );
        A.col(0) = Atmp.col(0)+0.5;
    }
}


////////////////////////////////////////////////////////////////////////////////
// Scaling the original aircraft coordinates to the grid. It is assumed that  //
// the input coods are not larger then one.                                   //
void CustomProfile::scalecoods(float scale, int gridsize, mat& A)
{
    A.col(0) = A.col(0)*scale + (gridsize-scale) /2;
    A.col(1) = A.col(1)*scale +  gridsize /2;
}

////////////////////////////////////////////////////////////////////////////////
// Creating the normal unit vectors on the original aircraft boundary
mat CustomProfile::createoriginalnormals (int whatcase, mat A, int Nrows)
{
    mat Normals = zeros(Nrows,2);
    int i = 0;
    mat N_magnitude = zeros(Nrows,1);

    if (whatcase == 1 )
    {
        while (i < Nrows)
        {
            if (i == 0)
            {
                Normals(i,0) =   ( A(i+1,1) - A(Nrows-1,1) );
                Normals(i,1) = - ( A(i+1,0) - A(Nrows-1,0) );
            }
            else if (i == Nrows-1)
            {
                Normals(i,0) =   ( A(0,1) - A(i-1,1) );
                Normals(i,1) = - ( A(0,0) - A(i-1,0) );
            }
            else
            {
                Normals(i,0) =   ( A(i+1,1) - A(i-1,1) );
                Normals(i,1) = - ( A(i+1,0) - A(i-1,0) );
            }
            i++;
        }
    }
    N_magnitude = sqrt ( Normals.col(0)%Normals.col(0) + Normals.col(1)%Normals.col(1));
    Normals.col(0) = Normals.col(0) / N_magnitude;
    Normals.col(1) = Normals.col(1) / N_magnitude;
    return(Normals);
}

////////////////////////////////////////////////////////////////////////////////
// Checking wether a point is inside a polygon.
int CustomProfile::pnpoly(int nvert, mat vertx, mat verty, float testx, float testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++)
  {
        if ( ((verty(i)>testy) != (verty(j)>testy)) &&	 (testx < (vertx(j)-vertx(i)) * (testy-verty(i)) / (verty(j)-verty(i)) + vertx(i)) )
        { c = !c; }
  }
  return c;
}

////////////////////////////////////////////////////////////////////////////////
// Creating the Exact Boundary Points, Normals and Arc lengths.
// Putting them in the format [x y ds nx ny]
mat CustomProfile::ExactBoundary( mat A, int Nrows)
{
    if ((A(Nrows-1,0) == A(1,0)) & (A(Nrows-1,1) == A(1,1))) Nrows--;
    mat ExactBndry = zeros( Nrows, 5);
    for (int i = 0;i<Nrows;i++)
    {
        if (i == Nrows -1 )    // For the last arc, calculate between last and first point.
        {
            if ((A(i,0) == A(1,0)) & (A(i,1) == A(1,1))) continue;  // If first and last point
                                                                // are the same: skip calculation
            ExactBndry(i,0) = (A(i,0) + A(0,0))/2;
            ExactBndry(i,1) = (A(i,1) + A(0,1))/2;
            ExactBndry(i,2) = sqrt( pow(A(0,0)-A(i,0),2.0) +
                                    pow(A(0,1)-A(i,1),2.0) );
            ExactBndry(i,3) =  (A(0,1)-A(i,1))/ExactBndry(i,2);
            ExactBndry(i,4) = - (A(0,0)-A(i,0))/ExactBndry(i,2);
            continue;
        }
        ExactBndry(i,0) = (A(i,0) + A(i+1,0))/2;
        ExactBndry(i,1) = (A(i,1) + A(i+1,1))/2;
        ExactBndry(i,2) = sqrt( pow(A(i+1,0)-A(i,0),2.0) +
                                pow(A(i+1,1)-A(i,1),2.0) );
        ExactBndry(i,3) =  (A(i+1,1)-A(i,1))/ExactBndry(i,2);
        ExactBndry(i,4) = - (A(i+1,0)-A(i,0))/ExactBndry(i,2);

    }
    return ExactBndry;
}

////////////////////////////////////////////////////////////////////////////////
// Creating the Grid Boundary Points and normals
// Putting them in the format [x y ds nx ny]
cube CustomProfile::CreateGridBndryCube(mat A, mat normals, int Nrows, int gridsizex, int gridsizey)
{
    mat Aint = round(A);
 //   mat NORMALSX = zeros(gridsizex,gridsizey), NORMALSY = zeros(gridsizex,gridsizey);
//    mat matGamma = zeros(gridsizex,gridsizey), Order = zeros(gridsizex,gridsizey);
    cube GridBndryCube(gridsizex,gridsizey,5);
//  Trying out sparse matrices. Couldnt write to file last time.
    sp_mat NORMALSX(gridsizex,gridsizey), NORMALSY(gridsizex,gridsizey);
    sp_mat matGamma(gridsizex,gridsizey), Order(gridsizex,gridsizey);
    int k=0;
    float j = 1;
    float i = 0, pointsbetweenX, pointsbetweenY, pointsbetween,x,y;
    while (i < Nrows-1)
    {
        pointsbetweenX = Aint(i+1,0)-Aint(i,0);
        pointsbetweenY = Aint(i+1,1)-Aint(i,1);
        pointsbetween  = max(abs(pointsbetweenX),abs(pointsbetweenY));
        if ((abs(pointsbetweenX)>= abs(pointsbetweenY)) & (pointsbetween >=1))
        {
            while (j <= pointsbetween)
            {
                matGamma(Aint(i,0)-1,Aint(i,1)-1)=1;
                x = round(   A(i,0)      +    j * pointsbetweenX/abs(pointsbetweenX) );
                y = round(  (A(i,1)*(pointsbetween-j) + A(i+1,1)*(j))/pointsbetween  );
                matGamma(x-1,y-1) = 1;
                Order(Aint(i,0)-1,Aint(i,1)-1)=k;
                Order(x-1,y-1) = k+j;
                NORMALSX(Aint(i,0)-1,Aint(i,1)-1)    =normals(i,0) +                                                           NORMALSX(Aint(i,0)-1,Aint(i,1)-1);
                NORMALSX(x-1,y-1)  =normals(i,0) *(pointsbetween-j)/pointsbetween + normals(i+1,0) * (j)/pointsbetween +       NORMALSX(x-1,y-1);
                NORMALSY(Aint(i,0)-1,Aint(i,1)-1)    =normals(i,1) +                                                           NORMALSY(Aint(i,0)-1,Aint(i,1)-1);
                NORMALSY(x-1,y-1)  =normals(i,1) *(pointsbetween-j)/pointsbetween + normals(i+1,1) * (j)/pointsbetween +       NORMALSY(x-1,y-1);
                j++;
            }
            k=k+j-1;
            j=1;
        }
        else if ((abs(pointsbetweenX)< abs(pointsbetweenY)) & (pointsbetween >=1))
        {
            while (j <= pointsbetween)
            {
                matGamma(Aint(i,0)-1,Aint(i,1)-1)=1;
                y = round(   A(i,1)      +    j * pointsbetweenY/abs(pointsbetweenY) );
                x = round(  (A(i,0)*(pointsbetween-j) + A(i+1,0)*(j))/pointsbetween  );
                matGamma(x-1,y-1) = 1;

                Order(Aint(i,0)-1,Aint(i,1)-1)=k;
                Order(x-1,y-1) = k+j;

                NORMALSX(Aint(i,0)-1,Aint(i,1)-1)    =normals(i,0) +                                                           NORMALSX(Aint(i,0)-1,Aint(i,1)-1);
                NORMALSX(x-1,y-1)  =normals(i,0) *(pointsbetween-j)/pointsbetween + normals(i+1,0) * (j)/pointsbetween +       NORMALSX(x-1,y-1);
                NORMALSY(Aint(i,0)-1,Aint(i,1)-1)    =normals(i,1) +                                                           NORMALSY(Aint(i,0)-1,Aint(i,1)-1);
                NORMALSY(x-1,y-1)  =normals(i,1) *(pointsbetween-j)/pointsbetween + normals(i+1,1) * (j)/pointsbetween +       NORMALSY(x-1,y-1);
                j++;
            }
            k=k+j-1;
            j=1;
        }
        else
        {

        }
        i++;
    }
    GridBndryCube.slice(0)=matGamma;
    GridBndryCube.slice(1)=Order;
    GridBndryCube.slice(2)=NORMALSX;
    GridBndryCube.slice(3)=NORMALSY;
    return GridBndryCube;
}

////////////////////////////////////////////////////////////////////////////////
// scaling down the matrix A to unit format. Such that: highest number found is 1
// It also sets the mean of A in y to half the screen.
// Be wise in when to call this function. It is not needed for airfoils or when the
// input format already corresponds the right gridsize.
mat CustomProfile::scalecoodsdown(mat A)
{
    float maxsize;
    mat maximums = max(A,0);
    mat minimums = min(A,0);
    if ((maximums(0)-minimums(0)) > (maximums(1)-minimums(1))) maxsize=(maximums(0)-minimums(0));
    else maxsize = (maximums(1)-minimums(1));
    A.col(0) = (A.col(0)-minimums(0))/(maximums(0)-minimums(0));
    A.col(1) = (A.col(1)-minimums(1)-maximums(1)/2)/(maximums(0)-minimums(0));
    return A;
}

////////////////////////////////////////////////////////////////////////////////
// Converting those big matrices to small arrays of Gamma [x,y] and NORMALS[Nx,Ny]
mat CustomProfile::CreateGridBndry(cube GridBndryCube, int nGamma, int gridsizex, int gridsizey)
{
    int k=0;
    mat matGamma=GridBndryCube.slice(0);
    mat Order = GridBndryCube.slice(1);
    mat NORMALSX=GridBndryCube.slice(2);
    mat NORMALSY=GridBndryCube.slice(3);
    mat gamma=zeros(nGamma,2);
    mat NORMALS = zeros(nGamma,2);
    vec ORDER = zeros(nGamma,1);
    for (int i=0; i<gridsizex; i++)              {
        for (int j=0; j<gridsizey; j++)          {
            if (matGamma(i,j)==1)               {
                gamma(k,0) = i;
                gamma(k,1) = j;
                ORDER(k) = Order(i,j);
                NORMALS(k,0) = NORMALSX(i,j);
                NORMALS(k,1) = NORMALSY(i,j);
                k++;            }               }
    }
    uvec indexorder=sort_index(ORDER);

    mat GridBndry = zeros(nGamma,4);
    for (int i = 0; i < nGamma; i++)
    {
        GridBndry(i,0)=gamma(indexorder(i),0);
        GridBndry(i,1)=gamma(indexorder(i),1);
        GridBndry(i,2) = NORMALS(indexorder(i),0);
        GridBndry(i,3) = NORMALS(indexorder(i),1);

    }
    // Set to unit vectors
    mat N_MAGNITUDE = zeros(nGamma,1);
    N_MAGNITUDE = sqrt (GridBndry.col(2) % GridBndry.col(2) + GridBndry.col(3) % GridBndry.col(3));
    GridBndry.col(2) = GridBndry.col(2)/N_MAGNITUDE;
    GridBndry.col(3) = GridBndry.col(3)/N_MAGNITUDE;
    return GridBndry;
}

////////////////////////////////////////////////////////////////////////////////
// Find interior boundary points Omega,using pnpoly from the web.
mat CustomProfile::FindOmega(mat GridBndry, int nGamma, int gridsizex, int gridsizey)
{
    int k=0, c;
//    mat omega=zeros(gridsizex,gridsizey);
    sp_mat omega(gridsizex,gridsizey);
    for (float i=0; i<gridsizex; i++)
    {
        for (float j=0; j<gridsizey; j++)
        {
            c = pnpoly(nGamma, GridBndry.col(0), GridBndry.col(1), i, j);
            if (c==1)
            {
                omega(i,j)=1;
            }
        }
    }

    // Remove boundary points from Omega
    for (float i=0; i<nGamma;i++) {omega(GridBndry(i,0),GridBndry(i,1))=0;}

    // Putting Omega in array format indicating x and y of each point
    int nOmega = accu(omega);
    mat Omega = zeros(nOmega,2);
    for (int i=0; i<gridsizex; i++)
    {
        for (int j=0; j<gridsizey; j++)
        {
            if (omega(i,j)==1)
            {
                Omega(k,0) = i;
                Omega(k,1) = j;
                k++;
            }
        }
    }
    return Omega;
}

void CustomProfile::setSize(vec x, vec y, float float1, float float2, float float3, float float4) {

}
