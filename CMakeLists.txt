#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TetCal)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#igl
#option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
#option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
#find_package(LIBIGL REQUIRED QUIET)

#Eigen
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#DQS-please add link to dqs folder
include_directories(${PROJECT_SOURCE_DIR}/dqs)
file(GLOB dqs ${PROJECT_SOURCE_DIR}/dqs/*)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(TetCal TetCal.cc ${sources} ${headers} ${dqs})
target_link_libraries(TetCal ${Geant4_LIBRARIES})
# Eigen3::Eigen igl::core)
# igl::opengl_glfw)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TetCal DESTINATION bin)
