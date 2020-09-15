cmake_minimum_required(VERSION 3.12)
project(iPref)

set(CMAKE_CXX_STANDARD 14)

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

add_executable(iPref collection.cpp
        collection.h
        filemem.cpp
        filemem.h
        global.cpp
        global.h
        header.h
        hypercube.cpp
        hypercube.h
        iPref.cpp
        iPref.h
        main.cpp
        mem.h
        param.cpp
        param.h
        point.cpp
        point.h
        queryGen.cpp
        queryGen.h
        rentry.cpp
        rentry.h
        rnode.cpp
        rnode.h
        rtree.cpp
        rtree.h
        skyline.cpp
        skyline.h
        tgs.cpp
        tgs.h
        virtualRNode.cpp
        virtualRNode.h
        qp_solver.h
        qp_solver_impl.h
        qhull_user.h
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
        lp_user.cpp
        lp_user.h)

target_link_libraries(${PROJECT_NAME} PRIVATE osqp::osqpstatic)
target_link_libraries(${PROJECT_NAME} PRIVATE libqhullcpp)

