#This CMakeLists is vaguely based off of Andrew Olivier's from https://github.com/MinervaExpt/NucCCNeutrons/blob/develop/CMakeLists.txt
#It's current intention is solely for the MINERvA GPVMs...

#CMake version set up by default on the MINERvA GPVMs
cmake_minimum_required(VERSION 2.8.12)

project(MINERvANeutronPurity)

#Copying Andrew's Compiler Flags (NOTE: STD of c++1y is worth keeping an eye on if things start behaving unexpectedly...)
#Also note that this defines the BUILD_TYPE it seems????????
set( GCC_Flags_For_CXX "-std=c++1y -Wall" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_Flags_For_CXX}" )
set( CMAKE_CXX_FLAGS_DEBUG "-ggdb" )
set( CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG" )

#Tell this package where it is installed and version control status
add_definitions(-DINSTALL_DIR="${CMAKE_INSTALL_PREFIX}/")

#Let directories in this package see each other
include_directories( "${PROJECT_SOURCE_DIR}" )

#Find dependencies
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree Core Geom EG GenVector OPTIONAL_COMPONENTS Cintex) #TODO: Remove CintEx if ROOT version > 5.34

if(${ROOT_VERSION} VERSION_LESS 6 AND NOT ${ROOT_Cintex_FOUND})
  MESSAGE(FATAL_ERROR "Cintex is optional except when it's not.  ROOT 6 has Reflex "
                      "support built in, so it doesn't have a separate component "
                      "for Cintex.  Cintex was an experimental feature in ROOT 5, "
                      "so I have to require it as a component there.  You appear to "
                      "be using ${ROOT_VERSION}, so I can't let you get away with "
                      "skipping out on Cintex support!  I need Cintex to load metadata "
                      "about MnvH1D and MnvH2D at runtime.")
endif()

if(${ROOT_VERSION} VERSION_LESS 6)
  MESSAGE("I need Cintex on this platform, so making sure to compile against it...")
else()
  MESSAGE("I don't need Cintex on this platform, so skipping it...")
  add_definitions(-DNCINTEX)
endif()

include(${ROOT_USE_FILE})

#Need to define the directory env variables $PlotUtils_DIR...
find_package(PlotUtils REQUIRED)
include_directories(${PlotUtils_INCLUDE_DIR})
message("Included PlotUtils from ${PlotUtils_INCLUDE_DIR}")

find_package(UnfoldUtils REQUIRED)
include_directories(${UnfoldUtils_INCLUDE_DIR})
message("Included UnfoldUtils from ${UnfoldUtils_INCLUDE_DIR}")

#Build main executables
add_executable(testExecutable testExecutable.cxx)

#link
target_link_libraries(testExecutable ${ROOT_LIBRARIES} PlotUtils UnfoldUtils)

#install
install(TARGETS testExecutable DESTINATION bin)
