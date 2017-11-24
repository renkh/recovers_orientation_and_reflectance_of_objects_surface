// Class for representing a 2D gray-scale image,
// with support for reading/writing pgm images.
// To be used in Computer Vision class.

#ifndef COMPUTER_VISION_IMAGE_H_
#define COMPUTER_VISION_IMAGE_H_

#include <cstdlib>
#include <string>
#include <fstream>

namespace ComputerVisionProjects {

// Class for representing a gray-scale image.
// Sample usage:
//   Image one_image;
//   one_image.AllocateSpaceAndSetSize(100, 200);
//   one_image.SetNumberGrayLevels(255);
//   // Creates and image such that each pixel is 150.
//   for (int i = 0; i < 100; ++i)
//     for (int j = 0; j < 200; ++j)
//       one_image.SetPixel(i, j, 150);
//   WriteImage("output_file.pgm", an_image);
//   // See image_demo.cc for read/write image.
class Image {
 public:
  Image(): num_rows_{0}, num_columns_{0},
	   num_gray_levels_{0}, pixels_{nullptr} { }

  Image(const Image &an_image);
  Image& operator=(const Image &an_image) = delete;

  ~Image();

  // Sets the size of the image to the given
  // height (num_rows) and columns (num_columns).
  void AllocateSpaceAndSetSize(size_t num_rows, size_t num_columns);

  size_t num_rows() const { return num_rows_; }
  size_t num_columns() const { return num_columns_; }
  size_t num_gray_levels() const { return num_gray_levels_; }
  void SetNumberGrayLevels(size_t gray_levels) {
    num_gray_levels_ = gray_levels;
  }

  size_t GetNumberGrayLevels() {
    return num_gray_levels_;
  }

  // Sets the pixel in the image at row i and column j
  // to a particular gray_level.
  void SetPixel(size_t i, size_t j, int gray_level) {
    if (i >= num_rows_ || j >= num_columns_) abort();
    pixels_[i][j] = gray_level;
  }

  int GetPixel(size_t i, size_t j) const {
    if (i >= num_rows_ || j >= num_columns_) abort();
    return pixels_[i][j];
  }

  int GetNumberOfRows() {
    return num_rows_;
  }

  int GetNumberOfColumns() {
    return num_columns_;
  }

 private:
  void DeallocateSpace();

  size_t num_rows_;
  size_t num_columns_;
  size_t num_gray_levels_;
  int **pixels_;
};

// Reads a pgm image from file input_filename.
// an_image is the resulting image.
// Returns true if  everyhing is OK, false otherwise.
bool ReadImage(const std::string &input_filename, Image *an_image);

// Writes image an_iamge into the pgm file output_filename.
// Returns true if  everyhing is OK, false otherwise.
bool WriteImage(const std::string &output_filename, const Image &an_image);

//  Draws a line of given gray-level color from (x0,y0) to (x1,y1);
//  an_image is the output_image.
// IMPORTANT: (x0,y0) and (x1,y1) can lie outside the image
//   boundaries, so SetPixel() should check the coordinates passed to it.
void DrawLine(int x0, int y0, int x1, int y1, int color,
	      Image *an_image);

/**
 * ConvertToBinary( ) sets image pixels to 0 if its value is below threshold
 * and 1 if its value is above threshold
 *
 * @param {int} threshold_value: the threshold value
 * @param {Image} an_image: input image
 */
void ConvertToBinary(const int threshold_value, Image *an_image);

/**
 * ComputeObjectAttributes( ) computes attributes that serve as object model
 * database. Atrributes include object label, row position of the center,
 * column position of the center, the minimum moment of inertia,
 * and the orientation
 *
 * @param {ostream} output_file: the output database file containing attributes
 * @param {Image} an_image: input image
 */
void ComputeObjectAttributes(std::ostream &output_file, Image *an_image);

/**
 * CompareObjectAttributes( ) compares the attributes of each object in a
 * labeled image file with those from the object model database.
 *
 * @param {istream} database_file: database file containing attributes
 * @param {Image} an_image: input image
 */
void CompareObjectAttributes(std::istream &database_file, Image *an_image);

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
void NormalToBrightestSurfaceSpot(const std::string input_file, std::ostream &output_file, Image *an_image);

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
void ComputeSurfaceNormals(const std::string input_file, Image *an_image_1, Image *an_image_2, Image *an_image_3, int step, int threshold);

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
 * @param threshold: threshold value for which only the pixels with intensity
 *                   greater than the threshold value will be counted.
 */
void ComputeSurfaceAlbedos(const std::string input_file, Image *an_image_1, Image *an_image_2, Image *an_image_3, int threshold);

void PrintImageToCout(Image *an_image);

}  // namespace ComputerVisionProjects

#endif  // COMPUTER_VISION_IMAGE_H_
