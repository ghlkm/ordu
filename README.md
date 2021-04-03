# iPref

Before running this code:

Step1:
download cmake

Step2:
download source code of Qhull, compile it using cmake and install it:http://www.qhull.org/download/

Step3:
modify the paths in CMakeLists.txt if your path of Qhull is different from us

To install qhull
-----------------
Installing Qhull with CMake 2.6 or later

  See CMakeLists.txt for examples and further build instructions

  To build Qhull, static libraries, shared library, and C++ interface
  - Download and extract Qhull (either GitHub, .tgz file, or .zip file)
  - cd build
  - cmake --help  # List build generators
  - make -G "<generator>" .. && cmake ..  
  - cmake ..
  - make
  - make install

  The ".." is important.  It refers to the parent directory (i.e., qhull/)

  On Windows, CMake installs to C:/Program Files/qhull.  64-bit generators
  have a "Win64" tag.  Qhull's data structures are substantial larger as
  64-bit code than as 32-bit code.  This may slow down Qhull.

  If creating a qhull package, please include a pkg-config file based on build/qhull*.pc.in

  If cmake fails with "No CMAKE_C_COMPILER could be found"
  - cmake was not able to find the build environment specified by -G "..."

-----------------



the source code for the bounded size top-k query (ORD, ORU) problem.

to run:

-w 1 -k 5 -d 4 -X 50 -f ./pdt/pdt4d1600k.txt -m UTK_BB -W ./user/user4d200k.txt -n 400000 -i ./idx/idx7.txt 