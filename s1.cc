/******************************************************************************
 * Title          : s1.cc
 * Author         : Renat Khalikov
 * Created on     : November 21, 2017
 * Description    : converts a gray–level image to a binary one using
 *                  a threshold value, and locates the sphere in an image and
 *                  computes its center and radius
 * Purpose        :
 * Usage          : ./s1 sphere0.pgm 100 parameters.txt
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
  if (argc!=4) {
    printf("Usage: %s {input original image} {input threshold value} {output parameters file}\n", argv[0]);
    return 0;
  }
  const string input_file(argv[1]);
  const string value(argv[2]);
  const string output_file(argv[3]);

  Image an_image;
  if (!ReadImage(input_file, &an_image)) {
    cout <<"Can't open file " << input_file << endl;
    return 0;
  }
  int threshold_value = stoi(value);  // convert string to int
  ConvertToBinary(threshold_value, &an_image);

  ofstream output_parameters(output_file);
  if (output_parameters.fail()) {
    cerr << "Could not open: {output database}\n";
    exit(1); // 1 indicates an error occurred
  }
  ComputeObjectAttributes(output_parameters, &an_image);
  output_parameters.close();
}
