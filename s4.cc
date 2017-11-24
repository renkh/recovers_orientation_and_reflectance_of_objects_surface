/******************************************************************************
 * Title          : s4.cc
 * Author         : Renat Khalikov
 * Created on     : November 21, 2017
 * Description    : computes the surface albedo
 * Purpose        :
 * Usage          : ./s4 directions.txt sphere1.pgm sphere2.pgm sphere3.pgm 100 sphere_s4_output.pgm
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
  if (argc!=7) {
    printf("Usage: %s {input directions} {image 1} {image 2} {image 3} {threshold} {output}\n", argv[0]);
    return 0;
  }
  const string input_file(argv[1]);
  const string image_1(argv[2]);
  const string image_2(argv[3]);
  const string image_3(argv[4]);
  const string threshold(argv[5]);
  const string output(argv[6]);

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

  ComputeSurfaceAlbedos(input_file, &an_image_1, &an_image_2, &an_image_3, stoi(threshold));

  if (!WriteImage(output, an_image_1)){
    cout << "Can't write to file " << output << endl;
    return 0;
  }
}
