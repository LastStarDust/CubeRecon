message("Cube Recon Library")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Configure the dependencies
find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# Define the source and include files that should be used for the library
# part of CubeRecon.
set(source
  CubeERepSim.cxx
  CubeInfo.cxx CubeHitUtilities.cxx CubeClusterManagement.cxx
  CubeMakeUsed.cxx CubeTimeSlice.cxx
  CubeHits3D.cxx CubeShareCharge.cxx
  CubeMakeHits3D.cxx
  )

set(includes
  CubeERepSim.hxx
  CubeInfo.hxx CubeHitUtilities.hxx CubeClusterManagement.cxx
  CubeMakeUsed.hxx CubeTimeSlice.hxx
  CubeHits3D.hxx CubeShareCharge.hxx
  CubeMakeHits3D.hxx
  )

# Make sure the current directories are available for the later
# compilation.  This is required to make ROOT_GENERATE_DICTIONARY
# work.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Build the library.
add_library(cuberecon SHARED ${source})

target_include_directories(cuberecon PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/CubeRecon>")

target_link_libraries(cuberecon PUBLIC cuberecon_io ${ROOT_LIBRARIES})

# Install the library for CubeRecon
install(TARGETS cuberecon
  EXPORT CubeReconTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/CubeRecon )

# Install the header files.
install(FILES ${includes} DESTINATION include/CubeRecon)
