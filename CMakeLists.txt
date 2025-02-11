#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(ecalib VERSION 1.0.0 LANGUAGES C CXX )

add_definitions(-std=c++20)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED COMPONENTS TMVA ROOTVecOps ROOTDataFrame Geom Eve)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(yaml-cpp REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${YAML_CPP_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include
		    ${ROOT_INCLUDE_DIRS})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ecalib main.cc ${sources} ${headers})
target_link_libraries(ecalib ${ROOT_LIBRARIES} yaml-cpp Minuit)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
# install(TARGETS ecalib DESTINATION bin)