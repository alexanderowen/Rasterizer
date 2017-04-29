#include <iostream>
#include <cmath>
#include <vtkDoubleArray.h>
#include <vtkDataSetWriter.h>
#include <vtkPointData.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

/************* 
 *
 * Math utility functions
 *
 *************/


/* Computes the cotangent */
double cot(double a)
{
   return 1. / tan(a); 
}

/* Traditional formula for linear interpolation. Returns nan if a == b */
double lerp(double a, double b, double fa, double fb, double x) 
{
    double t = (x - a) / (b - a);
    double l = fa + t * (fb - fa);
    return l;
}

/* Dot product between two 3D vectors */
double dotProduct(double a[], double b[]) 
{
    return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

/* Computes cross product between two 3D vectors and stores it in 'c' */
void crossProduct(double a[], double b[], double c[])
{
    c[0] = (a[1] * b[2]) - (a[2] * b[1]);
    c[1] = (b[0] * a[2]) - (a[0] * b[2]);
    c[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

/* Normalizes a 3D vector in place */
void normalize(double normal[]) 
{
    double length = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    if (length != 1.0) 
    {
        normal[0] = normal[0] / length;
        normal[1] = normal[1] / length;
        normal[2] = normal[2] / length;
    }
}

/* Custom double ceiling function; resolves approximation errors*/
double customCeil(double f)
{
    return ceil(f-0.00001);
}

/* Custom double floor function; resolves approximation errors*/
double customFloor(double f)
{
    return floor(f+0.00001);
}

/************* 
 *
 * VTK IO utility functions
 *
 *************/

vtkImageData *newImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void writeImage(vtkImageData *img, const char *filename)
{
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(filename);
    writer->Write();
    writer->Delete();
}


/************* 
 *
 * Phong shading model  
 *
 *************/

/* Specifies the lighting parameters applied to the view */
struct LightingParameters 
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] =  0.0;
         lightDir[2] = -0.8;
         Ka          =  0.3;
         Kd          =  0.7;
         Ks          =  5.3;
         alpha       =  7.5;
    };

    double lightDir[3];  // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};
LightingParameters lp;

/* Computes the shading factor as per the phong shading model */
double calculatePhongShading(LightingParameters &lp, double *viewDirection, double *normal) 
{
    double shadingFactor = 0.0;
    normalize(normal);
    normalize(lp.lightDir);
    // Ambient
    shadingFactor += lp.Ka;

    // Diffuse    
    shadingFactor += fabs(dotProduct(normal, lp.lightDir)) * lp.Kd;

    // Specular
    double h = 2.0 * (dotProduct(lp.lightDir, normal));
    double R[3];
    R[0] = (normal[0] * h) - lp.lightDir[0];
    R[1] = (normal[1] * h) - lp.lightDir[1];
    R[2] = (normal[2] * h) - lp.lightDir[2];

    double angle = dotProduct(viewDirection, R);
    double s     = pow(angle, lp.alpha);
    if (isnan(s)) // Could occur if lp.alpha is < 0
        shadingFactor += 0;
    else
        shadingFactor += lp.Ks * s;

    return shadingFactor;
} 

/* A Triangle is a basic graphical primitive; for rasterization, most operations will be
 * performed on it */
class Triangle
{
    public:
        double         X[3];
        double         Y[3];
        double         Z[3];
        double         colors[3][3];
        double         normals[3][3];
        double         shadingFactors[3];

        /* If no 2 vertices lie on the same y-line, triangle is irregular */
        int isIrregular() 
        {
            return Y[0] != Y[1] && Y[0] != Y[2] && Y[1] != Y[2];
        }

        /* If any 2 vertices are the same, triangle is collinear */
        int isCollinear() 
        {
            return (X[0] == X[1] && Y[0] == Y[1]) ||
                   (X[0] == X[2] && Y[0] == Y[2]) ||
                   (X[1] == X[2] && Y[1] == Y[2]); 
        }

        /* Return the index of the lowest vertex, based on y-value */
        int getMinVertexIndex() 
        {
            int index  = 0;
            for (int i = 1; i < 3; i++) 
                index = Y[i] < Y[index] ? i : index;
            return index;
        }

        /* Return the index of the highest vertex, based on y-value */
        int getMaxVertexIndex() 
        {
            int index  = 0;
            for (int i = 1; i < 3; i++) 
                index = Y[i] > Y[index] ? i : index; 
            return index;
        }
 
        /* Return the index of the 'middle' vertex index based on y-value */
        int getMiddleVertexIndex() 
        {
            int min = getMinVertexIndex(); 
            int max = getMaxVertexIndex();

            if ((min == 0 || min == 1) && (max == 0 || max == 1)) 
                return 2;
            else if ((min == 0 || min == 2) && (max == 0 || max == 2)) 
                return 1;
            else 
                return 0;
        }

        /* Normalize all normals of a triangle */
        void normalizeTriangle()
        {
            for (int i = 0; i < 3; i++)
                normalize(this->normals[i]);
        }

        /* Finds the intersection of the horizontal line at y = r and the line defined
         * by the ith and jth vertices. */ 
        double getXIntersection(int i, int j, double r) 
        {
            double slope = (Y[j] - Y[i]) / (X[j] - X[i]);
            if (std::isinf(slope))  // slope is vertical line 
                return X[i];        

            double b = Y[i] - (slope * X[i]);
            double intersect = (r - b) / slope;
            if (std::isnan(intersect))  // the two lines are the same 
                intersect =  X[i] < X[j] ? X[i] : X[j];
             
            return intersect;
        }

        /* Finds the intersection of the horizontal line at y = r
         * Determines the Z value by viewing the triangle from the z,y plane */
        double getZIntersection(int i, int j, double r) 
        {
            double slope = (Y[j] - Y[i]) / (Z[j] - Z[i]);
            if (std::isinf(slope))  // slope is vertical line 
                return Z[i];        

            double b = Y[i] - (slope * Z[i]);
            double intersect = (r - b) / slope;
            if (std::isnan(intersect))  // the two lines are the same 
                intersect =  Z[i] < Z[j] ? Z[i] : Z[j];
             
            return intersect;
        }
 
        /* Sets the horizontal line to lines[0]. Triangle should not be irregular. */
        void findHorizontal(int lines[3][2]) 
        {
            if (this->isIrregular()) {
                return; 
            }
            if (Y[0] == Y[1]) 
            {
                lines[0][0] = 0;
                lines[0][1] = 1;
                lines[1][0] = 0;
                lines[1][1] = 2;
                lines[2][0] = 1;
                lines[2][1] = 2;
            }           
            else if (Y[0] == Y[2]) 
            {
                lines[0][0] = 0;
                lines[0][1] = 2;
                lines[1][0] = 0;
                lines[1][1] = 1;
                lines[2][0] = 1;
                lines[2][1] = 2;
            }           
            else if (Y[1] == Y[2]) 
            {
                lines[0][0] = 1;
                lines[0][1] = 2;
                lines[1][0] = 0;
                lines[1][1] = 2;
                lines[2][0] = 0;
                lines[2][1] = 1;
            }           
        }

        /* Finds the lerp'd RGB values of x on the line 'line', and sets them to 'rgb' */
        void triangleColorLerp(double rgb[], int line[], double x, double altX) 
        {
            double a     = this->X[line[0]];
            double altA  = this->Y[line[0]];
            double b     = this->X[line[1]];
            double altB  = this->Y[line[1]];
            double *fa   = this->colors[line[0]];
            double *fb   = this->colors[line[1]];
            rgb[0] = lerp(a, b, fa[0], fb[0], x);
            if (isnan(rgb[0]))     // if a == b, this would occur
                rgb[0] = lerp(altA, altB, fa[0], fb[0], altX);
            rgb[1] = lerp(a, b, fa[1], fb[1], x);
            if (isnan(rgb[1])) 
                rgb[1] = lerp(altA, altB, fa[1], fb[1], altX);
            rgb[2] = lerp(a, b, fa[2], fb[2], x);
            if (isnan(rgb[2])) 
                rgb[2] = lerp(altA, altB, fa[2], fb[2], altX);
        }

        /* Finds the lerp'd normals of x on the line 'line, and sets them to 'normal' */
        void normalLerp(double normal[], int line[], double x, double altX) 
        {
            double a     = this->X[line[0]];
            double altA  = this->Y[line[0]];
            double b     = this->X[line[1]];
            double altB  = this->Y[line[1]];
            double *fa   = this->normals[line[0]];
            double *fb   = this->normals[line[1]];
            normal[0] = lerp(a, b, fa[0], fb[0], x);
            if (isnan(normal[0]))     // if a == b, this would occur
                normal[0] = lerp(altA, altB, fa[0], fb[0], altX);
            normal[1] = lerp(a, b, fa[1], fb[1], x);
            if (isnan(normal[1])) 
                normal[1] = lerp(altA, altB, fa[1], fb[1], altX);
            normal[2] = lerp(a, b, fa[2], fb[2], x);
            if (isnan(normal[2])) 
                normal[2] = lerp(altA, altB, fa[2], fb[2], altX);
        }

        /* Set the x,y,z values of dest */
        void setXYZandShading(Triangle *dest, int destI, int srcI) 
        {
            dest->X[destI] = this->X[srcI];
            dest->Y[destI] = this->Y[srcI];
            dest->Z[destI] = this->Z[srcI];
            dest->shadingFactors[destI] = this->shadingFactors[srcI];
        }

        /* Set the RGB values of dest */
        void setRGB(Triangle *dest, int destI, int srcI) 
        {
            dest->colors[destI][0] = this->colors[srcI][0];
            dest->colors[destI][1] = this->colors[srcI][1];
            dest->colors[destI][2] = this->colors[srcI][2];
        }

        /* Set the normal value of dest */
        void setNormals(Triangle *dest, int destI, int srcI) 
        {
            dest->normals[destI][0] = this->normals[srcI][0];
            dest->normals[destI][1] = this->normals[srcI][1];
            dest->normals[destI][2] = this->normals[srcI][2];
        }

        /* Split an irregular triangle into a flat-top triangle and flat-bottom triangle */
        void splitTriangle(Triangle *t1, Triangle *t2) 
        {
            Triangle t  = *this;
            int min = t.getMinVertexIndex();
            int max = t.getMaxVertexIndex();
            int mid = t.getMiddleVertexIndex();
        
            t.setXYZandShading(t1, 0, max);
            t.setXYZandShading(t1, 1, mid);
            t.setXYZandShading(t2, 0, min);
            t.setXYZandShading(t2, 1, mid);

            t.setRGB(t1, 0, max);
            t.setRGB(t1, 1, mid);
            t.setRGB(t2, 0, min);
            t.setRGB(t2, 1, mid);

            t.setNormals(t1, 0, max);
            t.setNormals(t1, 1, mid);
            t.setNormals(t2, 0, min);
            t.setNormals(t2, 1, mid);

            double newX = t.getXIntersection(max, min, t.Y[mid]);
            t1->X[2] = newX;
            t1->Y[2] = t.Y[mid];
            t2->X[2] = newX;
            t2->Y[2] = t.Y[mid];

            double newZ = getZIntersection(min, max, t.Y[mid]);
            t1->Z[2] = newZ;
            t2->Z[2] = newZ;

            double newRGB[3];
            int line[] = {min, max};
            t.triangleColorLerp(newRGB, line, t1->X[2], t1->Y[2]);
            t1->colors[2][0] = newRGB[0];
            t1->colors[2][1] = newRGB[1];
            t1->colors[2][2] = newRGB[2];
            t2->colors[2][0] = newRGB[0];
            t2->colors[2][1] = newRGB[1];
            t2->colors[2][2] = newRGB[2];

            double newNormals[3];
            t.normalLerp(newNormals, line, t1->X[2], t1->Y[2]);
            t1->normals[2][0] = newNormals[0];
            t1->normals[2][1] = newNormals[1];
            t1->normals[2][2] = newNormals[2];
            t2->normals[2][0] = newNormals[0];
            t2->normals[2][1] = newNormals[1];
            t2->normals[2][2] = newNormals[2];
            normalize(t1->normals[2]);
            normalize(t2->normals[2]);

            double newShadingFactor;
            newShadingFactor = lerp(X[min], X[max], shadingFactors[min], 
                    shadingFactors[max],t1->X[2] );
            if (isnan(newShadingFactor))
                newShadingFactor = lerp(Y[min], Y[max], shadingFactors[min], 
                        shadingFactors[max],t1->Y[2] );
            t1->shadingFactors[2] = newShadingFactor;
            t2->shadingFactors[2] = newShadingFactor;
        }
};

/* Reads a .vtk file and translates the data into triangle data  */
std::vector<Triangle> GetTriangles(char *filename)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName(filename);
    //cerr << "Reading" << endl;
    rdr->Update();
    //cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

/* Screen space class to which the mapping from geometry->image is performed  */
class Screen
{
    public:
        unsigned char   *colorBuffer;
        double          *depthBuffer;
        int              height;
        int              width;
        
        Screen(int w, int h, unsigned char *buf) 
        {
            width       = w;
            height      = h;
            colorBuffer = buf;
            depthBuffer = new double[height*width];
            for (int i = 0; i < height*width; i++) 
                depthBuffer[i] = -1.;
        }

       ~Screen() 
        {
            delete depthBuffer;
        }

        /* Initializes the buffers to default values */
        void initBuffers() 
        {
            for (int i = 0; i < height*width; i++) 
                depthBuffer[i] = -1.;
            for (int i = 0 ; i < height*width*3 ; i++)
                colorBuffer[i] = 0;
        }

        /* Sets the pixel color at coordinates (x,y) with the color (r,g,b) */ 
        void setColorBufferValue(int x, int y, double rgb[], double shadingValue) 
        {
            if (x >= width || y >= height || x < 0 || y < 0) 
                return;

            unsigned char *pixel = colorBuffer + (y * width * 3) + (x * 3);
            *(pixel)     = customCeil(fmin(1.0, shadingValue * rgb[0]) * 255);
            *(pixel + 1) = customCeil(fmin(1.0, shadingValue * rgb[1]) * 255);
            *(pixel + 2) = customCeil(fmin(1.0, shadingValue * rgb[2]) * 255);
        }

        /* Set value 'z' for the pixel (x, y) in the depth buffer*/
        void setDepthBufferValue(int x, int y, double z) 
        {
            if (x >= width || y >= height || x < 0 || y < 0) 
                return;
            depthBuffer[x + y*width] = z;
        }

        /* Return the depth buffer value at pixel (x, y) */
        double getDepthBufferValue(int x, int y) 
        {
            if (x >= width || y >= height || x < 0 || y < 0) 
                return -1;             // not a real z value
            return depthBuffer[x + y*width];
        }

        /* Perform the scanline algorithm using t */
        void scanlineAlgorithm(Triangle t) 
        {
            double rowMin = customCeil(t.Y[t.getMinVertexIndex()]);
            double rowMax = customFloor(t.Y[t.getMaxVertexIndex()]);
            for (int r = rowMin; r <= rowMax; r++) 
            {
                int lines[3][2];
                t.findHorizontal(lines);

                double intersect1 = t.getXIntersection(lines[1][0], lines[1][1], r);
                double intersect2 = t.getXIntersection(lines[2][0], lines[2][1], r);
                double leftX  = intersect1 < intersect2 ? intersect1 : intersect2;
                double rightX = intersect1 > intersect2 ? intersect1 : intersect2;

                int lLine[2], rLine[2];
                if (leftX == intersect1) 
                {
                    lLine[0] = lines[1][0];
                    lLine[1] = lines[1][1];
                    rLine[0] = lines[2][0];
                    rLine[1] = lines[2][1];
                } 
                else 
                {
                    lLine[0] = lines[2][0];
                    lLine[1] = lines[2][1];
                    rLine[0] = lines[1][0];
                    rLine[1] = lines[1][1];
                }
                double lRGB[3], rRGB[3];
                t.triangleColorLerp(lRGB, lLine, leftX, r);
                t.triangleColorLerp(rRGB, rLine, rightX, r);

                double lShadingFactor, rShadingFactor;
                lShadingFactor = lerp(t.X[lLine[0]], t.X[lLine[1]], t.shadingFactors[lLine[0]],
                        t.shadingFactors[lLine[1]], leftX);
                if (isnan(lShadingFactor))
                    lShadingFactor = lerp(t.Y[lLine[0]], t.Y[lLine[1]], t.shadingFactors[lLine[0]],
                        t.shadingFactors[lLine[1]], r);

                rShadingFactor = lerp(t.X[rLine[0]], t.X[rLine[1]], t.shadingFactors[rLine[0]],
                        t.shadingFactors[rLine[1]], rightX);
                if (isnan(rShadingFactor))
                    lShadingFactor = lerp(t.Y[rLine[0]], t.Y[rLine[1]], t.shadingFactors[rLine[0]],
                        t.shadingFactors[rLine[1]], r);

                double leftZ  = t.getZIntersection(lLine[0], lLine[1], r);
                double rightZ = t.getZIntersection(rLine[0], rLine[1], r);
                double z, cRGB[3];
                for (int c = customCeil(leftX); c <= customFloor(rightX); c++) 
                {
                    z           = lerp(leftX, rightX, leftZ, rightZ, c);
                    cRGB[0]     = lerp(leftX, rightX, lRGB[0], rRGB[0], c) ;
                    cRGB[1]     = lerp(leftX, rightX, lRGB[1], rRGB[1], c) ;
                    cRGB[2]     = lerp(leftX, rightX, lRGB[2], rRGB[2], c) ;
                    if (z > getDepthBufferValue(c, r)) 
                    {
                        double sf = lerp(leftX, rightX, lShadingFactor, rShadingFactor, c);
                        this->setColorBufferValue(c, r, cRGB, sf);
                        setDepthBufferValue(c, r, z);
                    }
                }
            }
        }
};

/* Class representing a 4x4 matrix */
class Matrix
{
  public:
        double A[4][4];

        /* Applies this matrix to 'ptIn' and saves result to 'ptOut'*/
        void transformPoint(const double *ptIn, double *ptOut)
        {
            ptOut[0] = ptIn[0]*A[0][0]
                     + ptIn[1]*A[1][0]
                     + ptIn[2]*A[2][0]
                     + ptIn[3]*A[3][0];
            ptOut[1] = ptIn[0]*A[0][1]
                     + ptIn[1]*A[1][1]
                     + ptIn[2]*A[2][1]
                     + ptIn[3]*A[3][1];
            ptOut[2] = ptIn[0]*A[0][2]
                     + ptIn[1]*A[1][2]
                     + ptIn[2]*A[2][2]
                     + ptIn[3]*A[3][2];
            ptOut[3] = ptIn[0]*A[0][3]
                     + ptIn[1]*A[1][3]
                     + ptIn[2]*A[2][3]
                     + ptIn[3]*A[3][3];
        }

        /* Given 2 matrices, compose them into 1 and return it */
        static Matrix composeMatrices(const Matrix &M1, const Matrix &M2)
        {
            Matrix rv;
            for (int i = 0 ; i < 4 ; i++)
                for (int j = 0 ; j < 4 ; j++)
                {
                    rv.A[i][j] = 0;
                    for (int k = 0 ; k < 4 ; k++)
                        rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
                }

            return rv;
        }
};



class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    /* Places the camera frame in 'c' in the form [v1, v2, v3, origin] */
    void getCameraFrame(double c[][3])
    {
        double OMF[3]; /* origin minus focus */
        OMF[0] = position[0] - focus[0];
        OMF[1] = position[1] - focus[1];
        OMF[2] = position[2] - focus[2];
        normalize(OMF);
        double v1[3], v2[3];
        
        crossProduct(up, OMF, v1);
        normalize(v1);
        crossProduct(OMF, v1, v2);
        normalize(v2);

        c[0][0] = v1[0];
        c[0][1] = v1[1];
        c[0][2] = v1[2];
        c[1][0] = v2[0];
        c[1][1] = v2[1];
        c[1][2] = v2[2];
        c[2][0] = OMF[0];
        c[2][1] = OMF[1];
        c[2][2] = OMF[2];
        c[3][0] = position[0];
        c[3][1] = position[1];
        c[3][2] = position[2];
    }

    Matrix viewTransform() 
    {
        Matrix m;
        m.A[0][0] = cot(angle / 2.);
        m.A[0][1] = 0.; 
        m.A[0][2] = 0.; 
        m.A[0][3] = 0.; 
        m.A[1][0] = 0.;
        m.A[1][1] = cot(angle / 2.);
        m.A[1][2] = 0.;
        m.A[1][3] = 0.;
        m.A[2][0] = 0.;
        m.A[2][1] = 0.;
        m.A[2][2] = (far + near)/(far - near);
        m.A[2][3] = -1.;
        m.A[3][0] = 0.;
        m.A[3][1] = 0.;
        m.A[3][2] = 2 * far * near / (far - near);
        m.A[3][3] = 0.;
        return m;
    }

    Matrix CameraTransform() 
    {
        double c[4][3];
        getCameraFrame(c); /* (v1, v2, v3, origin) */
        Matrix m;
        m.A[0][0] = c[0][0];
        m.A[0][1] = c[1][0];
        m.A[0][2] = c[2][0];
        m.A[0][3] = 0.;
        m.A[1][0] = c[0][1];
        m.A[1][1] = c[1][1];
        m.A[1][2] = c[2][1];
        m.A[1][3] = 0.;
        m.A[2][0] = c[0][2];
        m.A[2][1] = c[1][2];
        m.A[2][2] = c[2][2];
        m.A[2][3] = 0.;

        double t[3];
        t[0] = 0. - c[3][0];
        t[1] = 0. - c[3][1];
        t[2] = 0. - c[3][2];
        double tv[3], uv[3], ut[3];
        crossProduct(t, c[1], tv);
        crossProduct(c[0], c[1], uv);
        crossProduct(c[0], t, ut);
        m.A[3][0] = dotProduct(tv, c[2]) / dotProduct(uv, c[2]);
        m.A[3][1] = dotProduct(ut, c[2]) / dotProduct(uv, c[2]);
        m.A[3][2] = dotProduct(uv, t)    / dotProduct(uv, c[2]);
        m.A[3][3] = 1.;

        return m;
    }

    Matrix deviceTransform(double height, double width) 
    {
        Matrix m;
        m.A[0][0] = width / 2.;
        m.A[0][1] = 0.;
        m.A[0][2] = 0.;
        m.A[0][3] = 0.;
        m.A[1][0] = 0.;
        m.A[1][1] = height / 2.;
        m.A[1][2] = 0.;
        m.A[1][3] = 0.;
        m.A[2][0] = 0.;
        m.A[2][1] = 0.;
        m.A[2][2] = 1.;
        m.A[2][3] = 0.;
        m.A[3][0] = width / 2.;
        m.A[3][1] = height / 2.;
        m.A[3][2] = 0.;
        m.A[3][3] = 1.;
        return m;
    }
};

double sineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames - 2 * ramp;
    double height = 1. / (nNonRamp + 4 * ramp / M_PI);
    if (curFrame < ramp)
    {
        double factor = 2 * height * ramp / M_PI;
        double eval = cos(M_PI / 2 * ((double)curFrame) / ramp);
        return (1. - eval) * factor;
    }
    else if (curFrame > nFrames - ramp)
    {
        int amount_left = nFrames - curFrame;
        double factor = 2 * height * ramp / M_PI;
        double eval = cos(M_PI / 2 * ((double)amount_left / ramp));
        return 1. - (1 - eval) * factor;
    }
    double amount_in_quad = ((double)curFrame - ramp);
    double quad_part = amount_in_quad * height;
    double curve_part = height * (2 * ramp) / M_PI;
    return quad_part + curve_part;
}

Camera getCamera(int frame, int nframes)
{
    double t = sineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI / 6;
    c.position[0] = 40 * sin(2 * M_PI * t);
    c.position[1] = 40 * cos(2 * M_PI * t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

/* Calculates phong shading value for each Triangle in 'triangles' */
void calculateShadingFactors(std::vector<Triangle> *triangles, Camera c)   
{
    double view_dir[3];
    for (int i = 0; i < triangles->size(); i++) 
    {
        Triangle *t = &((*triangles)[i]);
        for (int j = 0; j < 3; j++) 
        {
            view_dir[0] = t->X[j] - c.position[0];
            view_dir[1] = t->Y[j] - c.position[1];
            view_dir[2] = t->Z[j] - c.position[2];
            normalize(view_dir);
            t->shadingFactors[j] = calculatePhongShading(lp, view_dir, t->normals[j]);
        }
    }
}

/* Transforms each Triangle in 't' from World space to Device space */
void transformTriangleToDeviceSpace(Triangle *t, Camera c, int height, int width)
{
    Matrix camTrans  = c.CameraTransform();
    Matrix viewTrans = c.viewTransform();
    Matrix devTrans  = c.deviceTransform(height, width); 
    Matrix totalTrans;
    totalTrans = Matrix::composeMatrices(camTrans, viewTrans);
    totalTrans = Matrix::composeMatrices(totalTrans, devTrans);
    for (int i = 0; i < 3; i++) 
    {
        double inPoints[] = {t->X[i], t->Y[i], t->Z[i], 1};
        double outPoints[4];
        totalTrans.transformPoint(inPoints, outPoints);
        t->X[i] = outPoints[0] / outPoints[3];
        t->Y[i] = outPoints[1] / outPoints[3];
        t->Z[i] = outPoints[2] / outPoints[3];
    }
}

int main(int argc, char* argv[])
{
    if (argc != 2) 
    {
        cerr << "Usage: " << argv[0] << " filename" << endl;
        return 1;
    }

    int width  = 1000;
    int height = 1000;
    vtkImageData *image   = newImage(width, height);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);

    std::vector<Triangle> tris = GetTriangles(argv[1]);
    Screen screen(width, height, buffer);

    int views[] = {0, 250, 500, 750};
    for (int i = 0; i < 4; i++) 
    {
        std::vector<Triangle> triangles = tris;
        screen.initBuffers();
        Camera c = getCamera(views[i], 1000);
        calculateShadingFactors(&triangles, c);

        for (int j = 0; j < triangles.size(); j++) 
        {
            Triangle t = triangles[j];
            transformTriangleToDeviceSpace(&t, c, height, width);

            if (t.isCollinear()) 
                continue;

            if (t.isIrregular()) 
            {
                Triangle t1, t2;
                t.splitTriangle(&t1,&t2);
                screen.scanlineAlgorithm(t1);
                screen.scanlineAlgorithm(t2);
            } 
            else 
                screen.scanlineAlgorithm(t);
        }
        char buf[256];
        sscanf(argv[1], "%*[^/]/%s", buf);
        char filename[256];
        sprintf(filename, "images/%s%d.png", buf, views[i]);
        fprintf(stdout, "Writing to %s\n", filename);
        writeImage(image, filename);
    }
    return 0;
}
