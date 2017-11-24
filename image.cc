
// Class for representing a 2D gray-scale image,
// with support for reading/writing pgm images.
// To be used in Computer Vision class.

#include "image.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>
#include <array>
#include <string>   // getline()

using namespace std;

namespace ComputerVisionProjects {

Image::Image(const Image &an_image){
  AllocateSpaceAndSetSize(an_image.num_rows(), an_image.num_columns());
  SetNumberGrayLevels(an_image.num_gray_levels());

  for (size_t i = 0; i < num_rows(); ++i)
    for (size_t j = 0; j < num_columns(); ++j){
      SetPixel(i,j, an_image.GetPixel(i,j));
    }
}

Image::~Image(){
  DeallocateSpace();
}

void
Image::AllocateSpaceAndSetSize(size_t num_rows, size_t num_columns) {
  if (pixels_ != nullptr) DeallocateSpace();
  pixels_ = new int*[num_rows];
  for (size_t i = 0; i < num_rows; ++i)
    pixels_[i] = new int[num_columns];

  num_rows_ = num_rows;
  num_columns_ = num_columns;
}

void
Image::DeallocateSpace() {
  for (size_t i = 0; i < num_rows_; i++)
    delete pixels_[i];
  delete pixels_;
  pixels_ = nullptr;
  num_rows_ = 0;
  num_columns_ = 0;
}

bool ReadImage(const string &filename, Image *an_image) {
  if (an_image == nullptr) abort();
  FILE *input = fopen(filename.c_str(),"rb");
  if (input == 0) {
    cout << "ReadImage: Cannot open file" << endl;
    return false;
  }

  // Check for the right "magic number".
  char line[1024];
  if (fread(line, 1, 3, input) != 3 || strncmp(line,"P5\n",3)) {
    fclose(input);
    cout << "ReadImage: Expected .pgm file" << endl;
    return false;
  }

  // Skip comments.
  do
    fgets(line, sizeof line, input);
  while(*line == '#');

  // Read the width and height.
  int num_columns,num_rows;
  sscanf(line,"%d %d\n", &num_columns, &num_rows);
  an_image->AllocateSpaceAndSetSize(num_rows, num_columns);


  // Read # of gray levels.
  fgets(line, sizeof line, input);
  int levels;
  sscanf(line,"%d\n", &levels);
  an_image->SetNumberGrayLevels(levels);

  // read pixel row by row.
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0;j < num_columns; ++j) {
      const int byte=fgetc(input);
      if (byte == EOF) {
        fclose(input);
        cout << "ReadImage: short file" << endl;
        return false;
      }
      an_image->SetPixel(i, j, byte);
    }
  }

  fclose(input);
  return true;
}

bool WriteImage(const string &filename, const Image &an_image) {
  FILE *output = fopen(filename.c_str(), "w");
  if (output == 0) {
    cout << "WriteImage: cannot open file" << endl;
    return false;
  }
  const int num_rows = an_image.num_rows();
  const int num_columns = an_image.num_columns();
  const int colors = an_image.num_gray_levels();

  // Write the header.
  fprintf(output, "P5\n"); // Magic number.
  fprintf(output, "#\n");  // Empty comment.
  fprintf(output, "%d %d\n%03d\n", num_columns, num_rows, colors);

  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_columns; ++j) {
      const int byte = an_image.GetPixel(i , j);
      if (fputc(byte,output) == EOF) {
	    fclose(output);
            cout << "WriteImage: could not write" << endl;
	    return false;
      }
    }
  }

  fclose(output);
  return true;
}

// Implements the Bresenham's incremental midpoint algorithm;
// (adapted from J.D.Foley, A. van Dam, S.K.Feiner, J.F.Hughes
// "Computer Graphics. Principles and practice",
// 2nd ed., 1990, section 3.2.2);
void
DrawLine(int x0, int y0, int x1, int y1, int color,
  Image *an_image) {
  if (an_image == nullptr) abort();

  #ifdef SWAP
  #undef SWAP
  #endif
  #define SWAP(a,b) {a^=b; b^=a; a^=b;}

  const int DIR_X = 0;
  const int DIR_Y = 1;

  // Increments: East, North-East, South, South-East, North.
  int incrE, incrNE, incrS, incrSE, incrN;
  int d;         /* the D */
  int x,y;       /* running coordinates */
  int mpCase;    /* midpoint algorithm's case */
  int done;      /* set to 1 when done */

  int xmin = x0;
  int xmax = x1;
  int ymin = y0;
  int ymax = y1;

  int dx = xmax - xmin;
  int dy = ymax - ymin;
  int dir;

  if (dx * dx > dy * dy) {  // Horizontal scan.
    dir=DIR_X;
    if (xmax < xmin) {
      SWAP(xmin, xmax);
      SWAP(ymin , ymax);
    }
    dx = xmax - xmin;
    dy = ymax - ymin;

    if (dy >= 0) {
      mpCase = 1;
      d = 2 * dy - dx;
    } else {
      mpCase = 2;
      d = 2 * dy + dx;
    }

    incrNE = 2 * (dy - dx);
    incrE = 2 * dy;
    incrSE = 2 * (dy + dx);
  } else {// vertical scan.
    dir = DIR_Y;
    if (ymax < ymin) {
      SWAP(xmin, xmax);
      SWAP(ymin, ymax);
    }
    dx = xmax - xmin;
    dy = ymax-ymin;

    if (dx >=0 ) {
      mpCase = 1;
      d = 2 * dx - dy;
    } else {
      mpCase = 2;
      d = 2 * dx + dy;
    }

    incrNE = 2 * (dx - dy);
    incrE = 2 * dx;
    incrSE = 2 * (dx + dy);
  }

  /// Start the scan.
  x = xmin;
  y = ymin;
  done = 0;

  while (!done) {
    an_image->SetPixel(x,y,color);

    // Move to the next point.
    switch(dir) {
      case DIR_X:
      if (x < xmax) {
       switch(mpCase) {
         case 1:
         if (d <= 0) {
          d += incrE;
          x++;
        } else {
          d += incrNE;
          x++;
          y++;
        }
        break;

        case 2:
        if (d <= 0) {
          d += incrSE;
          x++;
          y--;
        } else {
          d += incrE;
          x++;
        }
        break;
      }
    } else {
     done=1;
   }
   break;

   case DIR_Y:
   if (y < ymax) {
    switch(mpCase) {
     case 1:
     if (d <= 0) {
       d += incrE;
       y++;
     } else {
       d += incrNE;
       y++;
       x++;
     }
     break;

     case 2:
     if (d <= 0) {
      d += incrSE;
      y++;
      x--;
    } else {
      d += incrE;
      y++;
    }
    break;
	  } // mpCase
        } // y < ymin
        else {
         done=1;
       }
       break;
     }
   }
 }

/**
 * ConvertToBinary( ) sets image pixels to 0 if its value is below threshold
 * and 1 if its value is above threshold
 *
 * @param {int} threshold_value: the threshold value
 * @param {Image} an_image: input image
 */
void ConvertToBinary(const int threshold_value, Image *an_image) {
  if (an_image == nullptr) abort();
  int pixel;
  int row = an_image->GetNumberOfRows();
  int column = an_image->GetNumberOfColumns();
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < column; ++j) {
      pixel = an_image->GetPixel(i,j);
      if (pixel <= threshold_value)
        an_image->SetPixel(i,j,0);
      else
        an_image->SetPixel(i,j,1);
    }
  }
  an_image->SetNumberGrayLevels(1);
}

/**
 * ComputeObjectAttributes( ) computes attributes that serve as object model
 * database. Atrributes include object label, row position of the center,
 * column position of the center, the minimum moment of inertia,
 * and the orientation
 *
 * @param {ostream} output_file: the output database file containing attributes
 * @param {Image} an_image: input image
 */
void ComputeObjectAttributes(std::ostream &output_file, Image *an_image) {
  size_t number_of_obects = an_image->GetNumberGrayLevels();
  std::vector<int> area(number_of_obects+1, 0);
  std::vector<int> x_axis(number_of_obects+1, 0);
  std::vector<int> y_axis(number_of_obects+1, 0);
  std::vector<int> x_axis_squared(number_of_obects+1, 0);
  std::vector<int> y_axis_squared(number_of_obects+1, 0);
  std::vector<int> xy(number_of_obects+1, 0);

  int row = an_image->GetNumberOfRows();
  int column = an_image->GetNumberOfColumns();
  int leftmost = column;
  int rightmost = 0;
  int uppermost = row;
  int lowermost = 0;

  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < column; ++j) {
      int pixel = an_image->GetPixel(i,j);
      if (pixel != 0) {
        area[pixel] += 1;
        x_axis[pixel] += i;
        y_axis[pixel] += j;
        if (j < leftmost) leftmost = j;
        if (j > rightmost) rightmost = j;
        if (i < uppermost) uppermost = i;
        if (i > lowermost) lowermost = i;
      }
    }
  }
  // header
  output_file << "x-coordinate of the center | "
              << "y-coordinate of the center | "
              << "radius of the circle"
              << endl;

  for (int i = 1; i < number_of_obects+1; ++i) {
    // get the area and center of objects
    // center of the object (x, y)
    // x = (1/A)∑∑i bij
    // y = (1/A)∑∑j bij
    double y_pos_of_center = x_axis[i] / area[i];
    double x_pos_of_center = y_axis[i] / area[i];
    double diameter = ((rightmost - leftmost) + (lowermost - uppermost)) / 2;
    output_file << x_pos_of_center << " " << y_pos_of_center << " " << diameter/2 << endl;
  }
}

/**
 * CompareObjectAttributes( ) compares the attributes of each object in a
 * labeled image file with those from the object model database.
 *
 * @param {istream} database_file: database file containing attributes
 * @param {Image} an_image: input image
 */
void CompareObjectAttributes(istream &database_file, Image *an_image) {
  size_t number_of_obects = an_image->GetNumberGrayLevels();
  std::vector<int> area(number_of_obects+1, 0);
  std::vector<int> x_axis(number_of_obects+1, 0);
  std::vector<int> y_axis(number_of_obects+1, 0);
  std::vector<int> x_axis_squared(number_of_obects+1, 0);
  std::vector<int> y_axis_squared(number_of_obects+1, 0);
  std::vector<int> xy(number_of_obects+1, 0);

  int row = an_image->GetNumberOfRows();
  int column = an_image->GetNumberOfColumns();
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < column; ++j) {
      int pixel = an_image->GetPixel(i,j);
      if (pixel != 0) {
        area[pixel] += 1;
        x_axis[pixel] += i;
        y_axis[pixel] += j;
        x_axis_squared[pixel] += i*i;
        y_axis_squared[pixel] += j*j;
        xy[pixel] += i*j;
      }
    }
  }

  string header;
  getline(database_file, header);
  int database_object_label, database_row_center, database_column_center;
  double database_inertia, database_orientation;
  while(database_file.good()) {
    database_file >> database_object_label
                  >> database_row_center
                  >> database_column_center
                  >> database_inertia
                  >> database_orientation;

    for (int i = 1; i < number_of_obects+1; ++i) {
      // get the area and center of objects
      // center of the object (x, y)
      // x = (1/A)∑∑i bij
      // y = (1/A)∑∑j bij
      double x_pos_of_center = x_axis[i] / area[i];
      double y_pos_of_center = y_axis[i] / area[i];

      // calculate a, b and c using the formulas from the slides
      // a = ∫∫(x')^2 b(x, y)dx'dy'
      double a = x_axis_squared[i] - (area[i] * x_pos_of_center * x_pos_of_center);
      // 2∫∫(x'y') b(x, y)dx'dy'
      double b = 2 * (xy[i] - (area[i] * x_pos_of_center * y_pos_of_center));
      // c = ∫∫(y')2b(x, y)dx'dy'
      double c = y_axis_squared[i] - (area[i] * y_pos_of_center * y_pos_of_center);

      // Arc tangent of two numbers
      // The result is an angle expressed in radians.
      // To convert from radians to degrees, divide by 2
      // 2 x Pi radians = 360 degrees.
      double theta = atan2(b,a-c) / 2;
      // calculate E by using all the information obtained above
      // E = a sin^2(θ) − b sin(θ) cos(θ) + c cos^2(θ)
      double min_moment_of_inertia = a*sin(theta)*sin(theta) - b*sin(theta)*cos(theta) + c*cos(theta)*cos(theta);

      double smaller_number = min(min_moment_of_inertia, database_inertia);
      double larger_number = max(min_moment_of_inertia, database_inertia);
      double compare_inertia = smaller_number / larger_number;
      double threshold = 0.8;

      if (compare_inertia > 0.8) {
        int endpoint_x = x_pos_of_center + cos(theta)*50;
        int endpoint_y = y_pos_of_center + sin(theta)*50;
        DrawLine(x_pos_of_center, y_pos_of_center, endpoint_x, endpoint_y, 200, an_image);
      }
    }
  }
}

void PrintImageToCout(Image *an_image) {
  // matrix dimensions
  int row = an_image->GetNumberOfRows();
  int column = an_image->GetNumberOfColumns();
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < column; ++j) {
      std::cout << an_image->GetPixel(i,j) << " ";
    }
  }
}

/**
 * NormalToBrightestSurfaceSpot( ) compute the directions and intensities of
 * the light sources. Uses a formula, along with the parameters computed in
 * ./s1, to find the normal to the brightest surface spot on the sphere in each
 * of the 3 images. Assume that this is the direction of the corresponding
 * light source.
 *
 * @param input_file:  circle parameters, x and y coordinate of its center and
 *                     its radius.
 * @param output_file: directions file containing surface normal of an_image
 * @param an_image:    .pgm image used to compute its surface normal
 */
void NormalToBrightestSurfaceSpot(const std::string input_file, std::ostream &output_file, Image *an_image){
  ifstream parameters(input_file);
  if (parameters.fail()) {
    cerr << "Could not open: {input_file}\n";
    exit(1); // 1 indicates an error occurred
  }
  string line;
  getline(parameters, line); // kill header
  getline(parameters, line);
  stringstream number_stream(line);
  int x_center, y_center;
  double radius;
  number_stream >> x_center >> y_center >> radius;
  parameters.close();

  // estimate the highlight center
  // by finding the brightest pixel
  int brightest = -1;
  int x_brightest, y_brightest;
  int row = an_image->GetNumberOfRows();
  int column = an_image->GetNumberOfColumns();
  for (int i = 0; i < row; ++i){
    for (int j = 0; j < column; ++j){
      int pixel = an_image->GetPixel(i,j);
      if(brightest < pixel){
        brightest = pixel;
        x_brightest = j;
        y_brightest = i;
      }
    }
  }
  // use the formula to find (z-z0)
  // http://cosinekitty.com/raytrace/chapter06_sphere.html
  double x = x_brightest - x_center;
  double y = y_brightest - y_center;
  double z_value = sqrt(radius*radius - x*x - y*y);
  float normalize = - brightest / radius;
  output_file << x * normalize << " " << y * normalize << " " << z_value * normalize << endl;
}

/**
 * ComputeSurfaceNormals( ) given 3 images of an object, computes the normals
 * to that object’s surface. Pixel (x, y) is visible from all 3 light sources
 * if its brightness in all 3 images is greater than a certain threshold.
 * Compute the normals every N pixels along x and y axes. The value of N will
 * be supplied as step.
 *
 * @param input_file: directions file containing light source normals for
 *                    3 images with different light sources.
 * @param an_image_1: .pgm image with a particular light source
 * @param an_image_2: .pgm image with a particular light source
 * @param an_image_3: .pgm image with a particular light source
 * @param step:       value for which normals will be computed every "step"
 *                    pixels along the x and y axes.
 * @param threshold:  threshold value for which only the pixels with brightness
 *                    greater than the threshold value will be counted.
 */
void ComputeSurfaceNormals(const std::string input_file, Image *an_image_1, Image *an_image_2, Image *an_image_3, int step, int threshold){
  ifstream directions_file(input_file);
  if (directions_file.fail()) {
    cerr << "Could not open: {input_file}\n";
    exit(1); // 1 indicates an error occurred
  }

  // S-matrix
  std::vector<std::vector<double>> directions(3, std::vector<double>(3, 0));

  // input directions from input_file
  string line;
  for (int i = 0; i < 3; ++i){
    getline(directions_file, line);
    stringstream number_stream(line);
    for (int j = 0; j < 3; ++j){
      number_stream >> directions[i][j];
    }
  }

  // determinant needed for inverse formula
  double determinant = directions[0][0] * (directions[1][1] * directions[2][2] - directions[1][2] * directions[2][1]) -
                       directions[0][1] * (directions[1][0] * directions[2][2] - directions[1][2] * directions[2][0]) +
                       directions[0][2] * (directions[1][0] * directions[2][1] - directions[1][1] * directions[2][0]);

  // inverse S-matrix using inverse formula
  double inverse_x_image_1 = (directions[1][1] * directions[2][2] - directions[1][2] * directions[2][1]) / determinant;
  double inverse_y_image_1 = (directions[0][2] * directions[2][1] - directions[0][1] * directions[2][2]) / determinant;
  double inverse_z_image_1 = (directions[0][1] * directions[1][2] - directions[0][2] * directions[1][1]) / determinant;

  double inverse_x_image_2 = (directions[1][2] * directions[2][0] - directions[1][0] * directions[2][2]) / determinant;
  double inverse_y_image_2 = (directions[0][0] * directions[2][2] - directions[0][2] * directions[2][0]) / determinant;
  double inverse_z_image_2 = (directions[0][2] * directions[1][0] - directions[0][0] * directions[1][2]) / determinant;

  double inverse_x_image_3 = (directions[1][0] * directions[2][1] - directions[1][1] * directions[2][0]) / determinant;
  double inverse_y_image_3 = (directions[0][1] * directions[2][0] - directions[0][0] * directions[2][1]) / determinant;
  double inverse_z_image_3 = (directions[0][0] * directions[1][1] - directions[0][1] * directions[1][0]) / determinant;

  // assuming all images have same dimensions
  int row = an_image_1->GetNumberOfRows();
  int column = an_image_1->GetNumberOfColumns();
  for (int i = 0; i < row; i=i+step){
    for (int j = 0; j < column; j=j+step){
      int pixel_1 = an_image_1->GetPixel(i,j);
      int pixel_2 = an_image_2->GetPixel(i,j);
      int pixel_3 = an_image_3->GetPixel(i,j);
      if (pixel_1 >= threshold && pixel_2 >= threshold && pixel_3 >= threshold){
        // draw grid on an_image_1
        DrawLine(i, j, i, j, 0, an_image_1);

        // draw circle around each point
        DrawLine(i-1, j-1, i+1, j+1, 0, an_image_1);
        DrawLine(i-1, j+1, i+1, j-1, 0, an_image_1);

        // calculate surface normal to source normal
        double surface_normal_x = inverse_x_image_1 * pixel_1 + inverse_x_image_2 * pixel_2 + inverse_x_image_3 * pixel_3;
        double surface_normal_y = inverse_y_image_1 * pixel_1 + inverse_y_image_2 * pixel_2 + inverse_y_image_3 * pixel_3;
        double surface_normal_z = inverse_z_image_1 * pixel_1 + inverse_z_image_2 * pixel_2 + inverse_z_image_3 * pixel_3;
        double normalize = sqrt((surface_normal_x * surface_normal_x) + (surface_normal_y * surface_normal_y) + (surface_normal_z * surface_normal_z));

        // scale the vector 10 times after you normalize it to make it visible
        surface_normal_x = 10 * (surface_normal_x / normalize);
        surface_normal_y = 10 * (surface_normal_y / normalize);
        surface_normal_z = 10 * (surface_normal_z / normalize);

        DrawLine(i, j, i + surface_normal_x, j + surface_normal_y, 255, an_image_1);
      }
    }
  }
}

/**
 * ComputeSurfaceAlbedos( ) computes the surface albedo for all pixels visible
 * from all 3 light sources, scale them up or down to fit in the range 0...255
 * and show them in the output image.
 *
 * @param input_file: directions file containing light source normals for
 *                    3 images with different light sources.
 * @param an_image_1: .pgm image with a particular light source
 * @param an_image_2: .pgm image with a particular light source
 * @param an_image_3: .pgm image with a particular light source
 * @param threshold:  threshold value for which only the pixels with brightness
 *                    greater than the threshold value will be counted.
 */
void ComputeSurfaceAlbedos(const std::string input_file, Image *an_image_1, Image *an_image_2, Image *an_image_3, int threshold){
  ifstream directions_file(input_file);
  if (directions_file.fail()) {
    cerr << "Could not open: {input_file}\n";
    exit(1); // 1 indicates an error occurred
  }

  // S-matrix
  std::vector<std::vector<double>> directions(3, std::vector<double>(3, 0));

  // input directions from input_file
  string line;
  for (int i = 0; i < 3; ++i){
    getline(directions_file, line);
    stringstream number_stream(line);
    for (int j = 0; j < 3; ++j){
      number_stream >> directions[i][j];
    }
  }

  // determinant needed for inverse formula
  double determinant = directions[0][0] * (directions[1][1] * directions[2][2] - directions[1][2] * directions[2][1]) -
                       directions[0][1] * (directions[1][0] * directions[2][2] - directions[1][2] * directions[2][0]) +
                       directions[0][2] * (directions[1][0] * directions[2][1] - directions[1][1] * directions[2][0]);

  // inverse S-matrix using inverse formula
  double inverse_x_image_1 = (directions[1][1] * directions[2][2] - directions[1][2] * directions[2][1]) / determinant;
  double inverse_y_image_1 = (directions[0][2] * directions[2][1] - directions[0][1] * directions[2][2]) / determinant;
  double inverse_z_image_1 = (directions[0][1] * directions[1][2] - directions[0][2] * directions[1][1]) / determinant;

  double inverse_x_image_2 = (directions[1][2] * directions[2][0] - directions[1][0] * directions[2][2]) / determinant;
  double inverse_y_image_2 = (directions[0][0] * directions[2][2] - directions[0][2] * directions[2][0]) / determinant;
  double inverse_z_image_2 = (directions[0][2] * directions[1][0] - directions[0][0] * directions[1][2]) / determinant;

  double inverse_x_image_3 = (directions[1][0] * directions[2][1] - directions[1][1] * directions[2][0]) / determinant;
  double inverse_y_image_3 = (directions[0][1] * directions[2][0] - directions[0][0] * directions[2][1]) / determinant;
  double inverse_z_image_3 = (directions[0][0] * directions[1][1] - directions[0][1] * directions[1][0]) / determinant;

  // assuming all images have same dimensions
  int row = an_image_1->GetNumberOfRows();
  int column = an_image_1->GetNumberOfColumns();
  for (int i = 0; i < row; i++){
    for (int j = 0; j < column; j++){
      int pixel_1 = an_image_1->GetPixel(i,j);
      int pixel_2 = an_image_2->GetPixel(i,j);
      int pixel_3 = an_image_3->GetPixel(i,j);
      if (pixel_1 >= threshold && pixel_2 >= threshold && pixel_3 >= threshold){
        // calculate surface normal to source normal
        double surface_normal_x = inverse_x_image_1 * pixel_1 + inverse_x_image_2 * pixel_2 + inverse_x_image_3 * pixel_3;
        double surface_normal_y = inverse_y_image_1 * pixel_1 + inverse_y_image_2 * pixel_2 + inverse_y_image_3 * pixel_3;
        double surface_normal_z = inverse_z_image_1 * pixel_1 + inverse_z_image_2 * pixel_2 + inverse_z_image_3 * pixel_3;

        // calculate surface albedo
        double surface_albedo = sqrt((surface_normal_x * surface_normal_x) + (surface_normal_y * surface_normal_y) + (surface_normal_z * surface_normal_z));
        double surface_albedo_scaled = pixel_1 * surface_albedo;

       an_image_1->SetPixel(i,j,surface_albedo_scaled);
      }
    }
  }
}

}  // namespace ComputerVisionProjects
