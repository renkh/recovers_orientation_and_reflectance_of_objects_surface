/******************************************************************************
 * Title          : s2.cc
 * Author         : Renat Khalikov
 * Created on     : November 21, 2017
 * Description    : compute the directions and intensities of the light sources
 * Purpose        :
 * Usage          : ./s2 parameters.txt sphere1.pgm sphere2.pgm sphere3.pgm directions.txt
 * Build with     : make all
 */
#include "image.h"
#include <cstdio>
#include <iostream>
#include <fstream>  // ofstream
#include <string>

using namespace std;
using namespace ComputerVisionProjects;

int main(int argc, char **argv){
  if (argc!=6) {
    printf("Usage: %s {input parameters file} {image 1} {image 2} {image 3} {output directions file}\n", argv[0]);
    return 0;
  }
  const string input_file(argv[1]);
  const string image_1(argv[2]);
  const string image_2(argv[3]);
  const string image_3(argv[4]);
  const string output_file(argv[5]);

  Image an_image_1;
  if (!ReadImage(image_1, &an_image_1)) {
    cout <<"Can't open file " << image_1 << endl;
    return 0;
  }
  Image an_image_2;
  if (!ReadImage(image_2, &an_image_2)) {
    cout <<"Can't open file " << image_2 << endl;
    return 0;
  }
  Image an_image_3;
  if (!ReadImage(image_3, &an_image_3)) {
    cout <<"Can't open file " << image_3 << endl;
    return 0;
  }
  ofstream output_directions(output_file);
  if (output_directions.fail()) {
    cerr << "Could not open: {output directions file}\n";
    exit(1); // 1 indicates an error occurred
  }

  NormalToBrightestSurfaceSpot(input_file, output_directions, &an_image_1);
  NormalToBrightestSurfaceSpot(input_file, output_directions, &an_image_2);
  NormalToBrightestSurfaceSpot(input_file, output_directions, &an_image_3);

  output_directions.close();
}
