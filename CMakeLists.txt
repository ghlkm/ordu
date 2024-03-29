### the version of cmake can be lower or higher
cmake_minimum_required(VERSION 3.10)
project(iPref)

### the version of c++ can be later or equals to c++11
set(CMAKE_CXX_STANDARD 14)

### please see how to install qhull in README.md and change the path of qhull based on your machine
### add qhull
include_directories(/usr/local/include/libqhullcpp)

### qhullstatic reentrant library
add_library(libqhullstatic_r STATIC IMPORTED)
set_property(TARGET libqhullstatic_r PROPERTY IMPORTED_LOCATION /usr/local/lib/${CMAKE_STATIC_LIBRARY_PREFIX}qhullstatic_r${CMAKE_STATIC_LIBRARY_SUFFIX})
add_dependencies(libqhullstatic_r qhull)
### qhullcpp library
add_library(libqhullcpp STATIC IMPORTED)
set_property(TARGET libqhullcpp PROPERTY IMPORTED_LOCATION /usr/local/lib/${CMAKE_STATIC_LIBRARY_PREFIX}qhullcpp${CMAKE_STATIC_LIBRARY_SUFFIX})
set_property(TARGET libqhullcpp PROPERTY INTERFACE_LINK_LIBRARIES libqhullstatic_r)
add_dependencies(libqhullcpp qhull)

### please see how to install osqp in README.md and change the path of qhull based on your machine
### add osqp
find_package(osqp REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed -O3 -Wall -ldl")

include_directories(${CMAKE_SOURCE_DIR}/lp_lib/bfp)
include_directories(${CMAKE_SOURCE_DIR}/lp_lib/bfp/bfp_LUSOL)
include_directories(${CMAKE_SOURCE_DIR}/lp_lib/bfp/bfp_LUSOL/LUSOL)
include_directories(${CMAKE_SOURCE_DIR}/lp_lib/colamd)
include_directories(${CMAKE_SOURCE_DIR}/lp_lib/lp_solve)
include_directories(${CMAKE_SOURCE_DIR}/lp_lib/lpsolve55)
include_directories(${CMAKE_SOURCE_DIR}/lp_lib/shared)
include_directories(${CMAKE_SOURCE_DIR}/lp_lib)

include_directories(.)

add_executable(iPref
        iPref.cpp
        iPref.h
        main.cpp
        skyline.cpp
        skyline.h
        qp_solver.h
        qp_solver_impl.h
        qhull_user.h
        lp_user.cpp
        lp_user.h
        qp_solver2.cpp
        qp_solver.h
        case_study.cpp
        case_study.h

        rtree_lib/collection.cpp
        rtree_lib/collection.h
        rtree_lib/filemem.cpp
        rtree_lib/filemem.h
        rtree_lib/global.cpp
        rtree_lib/global.h
        rtree_lib/header.h
        rtree_lib/hypercube.cpp
        rtree_lib/hypercube.h
        rtree_lib/mem.h
        rtree_lib/param.cpp
        rtree_lib/param.h
        rtree_lib/point.cpp
        rtree_lib/point.h
        rtree_lib/rentry.cpp
        rtree_lib/rentry.h
        rtree_lib/rnode.cpp
        rtree_lib/rnode.h
        rtree_lib/rtree.cpp
        rtree_lib/rtree.h
        rtree_lib/tgs.cpp
        rtree_lib/tgs.h
        rtree_lib/virtualRNode.cpp
        rtree_lib/virtualRNode.h

        lp_lib/bfp/bfp_LUSOL/lp_LUSOL.h
        lp_lib/bfp/bfp_LUSOL/LUSOL/lusol.h
        lp_lib/bfp/lp_BFP.h
        lp_lib/colamd/colamd.h
        lp_lib/ini.h
        lp_lib/shared/commonlib.h
        lp_lib/shared/mmio.h
        lp_lib/shared/myblas.h
        lp_lib/lp_crash.h
        lp_lib/lp_Hash.h
        lp_lib/lp_lib.h
        lp_lib/lp_matrix.h
        lp_lib/lp_MDO.h
        lp_lib/lp_mipbb.h
        lp_lib/lp_MPS.h
        lp_lib/lp_presolve.h
        lp_lib/lp_price.h
        lp_lib/lp_pricePSE.h
        lp_lib/lp_report.h
        lp_lib/lp_scale.h
        lp_lib/lp_simplex.h
        lp_lib/lp_SOS.h
        lp_lib/lp_utils.h
        lp_lib/lp_wlp.h
        lp_lib/bfp/bfp_LUSOL/lp_LUSOL.c
        lp_lib/bfp/bfp_LUSOL/LUSOL/lusol.c
        lp_lib/colamd/colamd.c
        lp_lib/ini.c
        lp_lib/shared/commonlib.c
        lp_lib/shared/mmio.c
        lp_lib/shared/myblas.c
        lp_lib/lp_crash.c
        lp_lib/lp_Hash.c
        lp_lib/lp_lib.c
        lp_lib/lp_matrix.c
        lp_lib/lp_MDO.c
        lp_lib/lp_mipbb.c
        lp_lib/lp_MPS.c
        lp_lib/lp_params.c
        lp_lib/lp_presolve.c
        lp_lib/lp_price.c
        lp_lib/lp_pricePSE.c
        lp_lib/lp_report.c
        lp_lib/lp_scale.c
        lp_lib/lp_simplex.c
        lp_lib/lp_SOS.c
        lp_lib/lp_utils.c
        lp_lib/lp_wlp.c

        quadProg_lib/Array.hh
        quadProg_lib/Array.cc
        quadProg_lib/QuadProg++.hh
        quadProg_lib/QuadProg++.cc
        region.h)

target_link_libraries(${PROJECT_NAME} PRIVATE osqp::osqpstatic)
target_link_libraries(${PROJECT_NAME} PRIVATE libqhullcpp)

