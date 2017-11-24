Assignment
----------------
All parts of this assignment are complete.
Formula for s1 is z-z0 = sqrt(r^2 - (x-x0)^2 + (y-y0)^2)
Threshold value for s1 is 100
Step value for s3 is 10 and threshold value is 100
Threshold value for s4 is 100
----------------
Bugs:
Possible bug for albedo map, it comes out weird looking.
----------------
Compilation Instructions
----------------
Open a terminal.
Go to the directory containing the source code.
To compile type:
  make all
----------------
It is assumed that you are using a Linux machine and a g++ compiler.
----------------
To execute:
For s1
./s1 sphere0.pgm 100 parameters.txt

For s2
./s2 parameters.txt sphere1.pgm sphere2.pgm sphere3.pgm directions.txt

For s3
./s3 directions.txt sphere1.pgm sphere2.pgm sphere3.pgm 10 100 sphere_s3_output.pgm

For s4
./s4 directions.txt sphere1.pgm sphere2.pgm sphere3.pgm 100 sphere_s4_output.pgm
---------------
Note:
Prob. 1Libertian intensity is not dependent on the direction of view point, thus the ambient and diffuse intensities are independent as well. So the view point will not change the illumination of the sphere effectively causing the illumination to come from one direction s3. S3 is related to s1 and s2, as the intensity of s1 or s2 changes, diffuse intensities and ambient light will change as well causing s3 to change with it. If two distant light sources have unequal intensities, the direction of intensity of the effective source will point towards the light source with the higher intensity because the ambient light will be higher in that direction. 