# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
# File: CMakeLists.txt                                        #
#                                                             #
# - CMake configuration file to compile GiNaCDE V1.6.0          #
#                                                             #
# Project name: GiNaCDE - GiNaC Differential Equation Solver  #
# Contact: Mithun Bairagi <bairagirasulpur@gmail.com>         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 1.5.0                                                #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


cmake_minimum_required(VERSION 3.1)
project(GINACDE)

INCLUDE(GNUInstallDirs)

# information about the project
set( PROJECT_FULLNAME "GiNaCDE")
set( PROJECT_DESCRIPTION "GiNaCDE- an NLPDE solver")

# NOTE: On Windows and Cygwin, the dll's are placed in the
# CMAKE_RUNTIME_OUTPUT_DIRECTORY, while their corresponding import
# libraries end up in CMAKE_ARCHIVE_OUTPUT_DIRECTORY. On other
# platforms, such as Linux, the shared libraries are put in
# CMAKE_ARCHIVE_OUTPUT_DIRECTORY instead.
# Static libraries end up in CMAKE_ARCHIVE_OUTPUT_DIRECTORY on all
# platforms.
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)



# The version number.
set (GINACDE_VERSION_MAJOR 1)
set (GINACDE_VERSION_MINOR 6)
set(GINACDE_VERSION_MAINTENANCE 0)

set( GINACDE_VERSION "${GINACDE_VERSION_MAJOR}.${GINACDE_VERSION_MINOR}.${GINACDE_VERSION_MAINTENANCE}")
set( PROJECT_VERSION ${GINACDE_VERSION} )
set( GINACDE_NAME \"${PROJECT_FULLNAME}\")


# configure file specification in order to translate CMake variables to C++
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
set( DOLLAR "$" )
configure_file( ${CMAKE_SOURCE_DIR}/version.h.in ${CMAKE_BINARY_DIR}/version.h )
configure_file( ${CMAKE_SOURCE_DIR}/GiNaCDE.pc.in ${CMAKE_BINARY_DIR}/GiNaCDE.pc )
include_directories( ${CMAKE_BINARY_DIR} )

# specify pkg-config file to be installed
install( FILES ${CMAKE_BINARY_DIR}/GiNaCDE.pc DESTINATION lib/pkgconfig )

# add the binary tree to the search path for include files
# so that we will find version.h
include_directories("${PROJECT_BINARY_DIR}")

OPTION(GINACDE_BUILD_SHARED_LIBS "Build GiNaCDE shared library." ON)
OPTION(GINACDE_BUILD_STATIC_LIBS "Build GiNaCDE static library." ON)
OPTION(GINACDE_GUI_BUILD "Build GiNaCDE gui." OFF)

set (CMAKE_CXX_STANDARD 11)

# path to find own modules
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

#for finding external libraries ginac and cln.
find_package( GiNaC 1.7.6 REQUIRED )
find_package( cln 1.3.4 REQUIRED )

set(CORELIBS ${CLN_LIBRARIES} ${GINAC_LIBRARY})


include_directories(${GINACDE_SOURCE_DIR}/src)
link_directories(${GINACDE_BINARY_DIR}/lib)

SET(GINACDE_SRCS
    src/derivative.cpp
    src/desolve.cpp
    src/F_expns_methd.cpp
    src/fim.cpp
    src/inifcns.cpp
    src/integrate.cpp
    src/outform.cpp
    src/simplify.cpp
    src/solve.cpp
    src/utility.cpp
)
SET(GINACDE_HEADERS
    src/derivative.h
    src/desolve.h
    src/F_expns_methd.h
    src/fim.h
    src/inifcns.h
    src/integrate.h
    src/outform.h
    src/simplify.h
    src/solve.h
    src/utility.h
    src/GiNaCDE.h
)


# Building shared or static library
IF(GINACDE_BUILD_SHARED_LIBS)
    ADD_LIBRARY(GiNaCDE SHARED ${GINACDE_SRCS} ${GINACDE_HEADERS})
    target_link_libraries(GiNaCDE ${CORELIBS})
ENDIF()
IF(GINACDE_BUILD_STATIC_LIBS)
    ADD_LIBRARY(GiNaCDE_static STATIC ${GINACDE_SRCS} ${GINACDE_HEADERS})
    target_link_libraries(GiNaCDE_static ${CORELIBS})
ENDIF()

# Installing library
IF(GINACDE_BUILD_SHARED_LIBS)
    install (TARGETS GiNaCDE
             RUNTIME DESTINATION bin 
             LIBRARY DESTINATION lib
	     ARCHIVE DESTINATION lib)
ENDIF()
IF(GINACDE_BUILD_STATIC_LIBS)
    install (TARGETS GiNaCDE_static DESTINATION lib)
ENDIF()

# Testing library
include(CTest)


IF(GINACDE_BUILD_SHARED_LIBS)
    set(testlibs ${CORELIBS} GiNaCDE)
ENDIF()
IF(GINACDE_BUILD_STATIC_LIBS)
    set(testlibs ${CORELIBS} GiNaCDE_static)
ENDIF()

add_subdirectory(examples EXCLUDE_FROM_ALL)

add_executable(test1 test/test1.cpp)
target_link_libraries(test1 ${testlibs})

add_executable(test2 test/test2.cpp)
target_link_libraries(test2 ${testlibs})

add_executable(checkSolutest1 test/checkSolu_cahnAllen_mF.cpp)
target_link_libraries(checkSolutest1 ${testlibs})

add_executable(checkSolutest2 test/checkSolu_Painlev_FIMextravar.cpp)
target_link_libraries(checkSolutest2 ${testlibs})

add_test(test1Run ${GINACDE_BINARY_DIR}/bin/test1)
add_test(test2Run ${GINACDE_BINARY_DIR}/bin/test2)
add_test(checkSolutest1Run ${GINACDE_BINARY_DIR}/bin/checkSolutest1)
add_test(checkSolutest2Run ${GINACDE_BINARY_DIR}/bin/checkSolutest2)

# Building gui program
add_executable(gtools tools/gtools.cpp)
target_link_libraries(gtools ${CORELIBS} GiNaCDE)

IF(GINACDE_GUI_BUILD)
    find_package (PkgConfig REQUIRED)
    pkg_check_modules (GTK3 REQUIRED gtk+-3.0)
    include_directories (${GTK3_INCLUDE_DIRS})
    link_directories (${GTK3_LIBRARY_DIRS})
    add_executable(GiNaCDE_gui tools/GiNaCDE_gui.cpp src/icon.rc ${GINACDE_SRCS} ${GINACDE_HEADERS})
    target_compile_definitions(GiNaCDE_gui PUBLIC -DGiNaCDE_gui)
    target_link_libraries(GiNaCDE_gui ${GTK3_LIBRARIES} ${CORELIBS})
ENDIF()

install (TARGETS gtools DESTINATION bin)

IF(GINACDE_GUI_BUILD)
install (TARGETS GiNaCDE_gui DESTINATION bin)
ENDIF()
# Installing headers
install (FILES ${GINACDE_HEADERS} DESTINATION include/GiNaCDE)
install (FILES "${PROJECT_BINARY_DIR}/version.h" DESTINATION include/GiNaCDE)
