# iPref

WE HAVE UPLOAD THE BASIC NON-ORDER SENSITIVE ORU OPERATOR.

The code of this project is cleaning

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
  - git clone https://github.com/qhull/qhull.git
  - cd qhull/build
  - cmake --help  # List build generators
  - cmake -G "<generator>" ..   
  - cmake --build .
  - sudo cmake --build . --target install


  The ".." is important.  It refers to the parent directory (i.e., qhull/)

  On Windows, CMake installs to C:/Program Files/qhull.  64-bit generators
  have a "Win64" tag.  Qhull's data structures are substantial larger as
  64-bit code than as 32-bit code.  This may slow down Qhull.

  If creating a qhull package, please include a pkg-config file based on build/qhull*.pc.in

  If cmake fails with "No CMAKE_C_COMPILER could be found"
  - cmake was not able to find the build environment specified by -G "..."

-----------------

To install osqp
-----------------

  - git clone --recursive https://github.com/oxfordcontrol/osqp
  - cd osqp
  - mkdir build
  - cd build
  - cmake -G "Unix Makefiles" .. #
  - cmake --build .
  - sudo cmake --build . --target install  
-----------------

  
  
the source code for the bounded size top-k query (ORD, ORU) problem.

to run:

-w 10 -k 10 -d 6 -m 90 -f ./pdt/HOUSE6D.dat -mt CS5 -W ./user/user6d200k.txt -n 2000000 -i ./idx/idx4.txt
