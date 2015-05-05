
cmake_minimum_required(VERSION 2.8)

project(MRI_RECON)
set(PROGNAME recon_binary)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CmakeModules/")
file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

ADD_DEFINITIONS(-DLINUX -DRECON_64BIT -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -Wall -Wextra -m64 -J16 -O2 -fopenmp")
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} )

# Find Required Libraries
message( STATUS "Check for required libraries ")
include(FindBlitz)
include(FindFFTW)
include(FindHDF5)
include(FindARMADILLO)

#Set include directories
include_directories( BEFORE ${ARMADILLO_INCLUDE_DIR}  ${BLITZ_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} )
set(LOCAL_LIBS ${FFTW3_LIBRARIES} ${ARMADILLO_LIBRARIES} ${BLITZ_LIBRARIES} ${HDF5_LIBRARIES} dl m z )

set(CMAKE_FIND_STATIC FIRST) 
set(BUILD_SHARED_LIBS ON)
set(STATIC_LIBS ON)

link_directories(${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH} )


#--------------------------------------
# Support for Command Line Binary
#--------------------------------------

add_library(SUPPORTLIB STATIC gridFFT.cpp ArrayTemplates.cpp )
target_link_libraries(SUPPORTLIB ${LOCAL_LIBS} )


add_executable(recon_binary recon.cxx ${LOCAL_INCS})
target_link_libraries(recon_binary SUPPORTLIB ${LOCAL_LIBS} )










