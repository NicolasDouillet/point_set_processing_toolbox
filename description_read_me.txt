Please first check the tests section get some relevant and use examples of this toobox functions.

To use the sources, mind to add its path to your paths, just like in the header of any test file : addpath(genpath('my_path_to_the_sources'));

Please don'tforget ta rate if this code helped you. Thanks ! ;-)


%% DESCRIPTION

This package is a point set processing toolbox which uses command lines in Matlab(R) console. It is designed to deal with and process 3D point sets with or without colors.
For 2D point set processing, you can just padd your P = [X,Y] coordinates with one Z zeros vector : P = cat(2,P,zeros(size(P,1),1)]


%% HELP

A basic help is included in the header of each source file. It especially includes the function role and its input and output arguments precise descriptions (role, class, size, etc).
Just like for any Matlab (R) function, typewrite "help my_mesh_processing_file" in Matlab console to access it.


%% DATA FORMATS & HYPOTHESIS ON THE POINT SET

Most of the functions included take very common and widely used data structures as inputs and outputs :

- P : the point set / point cloud. Real matrix of double. size(P) = [nb_pts,3].
- N : the point set normals.       Real matrix of double. size(N) = [nb_pts,3].
- C : the color set. 	           Positiv real matrix of double. size(C) = [nb_pts,3].

where :

- A point is a 1 x 3 row double vector of real numbers.

By default, and unless exceptions, point arguments are index based.

Another common argument is nb_ngb which corresponds to the number of neighbors whished and used for the choosen process.
For commun usages, this value is in the range |[4 ; ~50]|. Tune it relatively to the local density of your point set


%% TESTS

Run the test files to discover and test the possibilities of this toolbox.
Use .mat data files provided in /data for test and example files.
Most of the functions and every important ones have been tested in a dedicated file named : test_my_function.m.
Note that no mesh reader or writer is provided in this toolbox since there already exist enough satisfying ones coded in Matlab.
Look for : read_ply.m, write_ply.m, read_off.m, write_off.m, plyread.m, plywrite.m
Then to create your own .mat file for point set P, just use the command save('path_to_my_file/my_file.mat','P');
To use the sources, mind to add its path to your paths, just like in the header of any test file : addpath(genpath('my_path_to_the_sources'));


%% COPYRIGHT, LICENSE AND TERMS OF USE

Please read the attached file : Copyright_license_&_terms_of_use.txt


%% MISC INFORMATION

By default, point normals are normalized and oriented at the same time they are computed.

Basic 3D mathematical computation algorithm (like point_to_plane_distance in src/geometry_ressources) are also independently available with their documentations in my file exchange contributions.

Matlab users, your advices and tips to improve and speed up my algorithms are welcome !

Since I am not native english speaker, please forgive my langage approximations.

Matlab release version used for development and tests : R2019b.

Contact : please report me any bug (with data set used and Matlab(R) code attached) or suggestion at nicolas.douillet9 (at) gmail.com

Last update : 27 / 12 /2024.