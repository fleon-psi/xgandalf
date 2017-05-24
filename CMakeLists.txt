cmake_minimum_required(VERSION 3.1)
project(indexer)
 
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo")
endif()
 
include_directories(include)

set(SOURCES src/Dbscan.cpp
            src/DetectorToReciprocalSpaceTransform.cpp
            src/ExperimentSettings.cpp
            src/HillClimbingOptimizer.cpp
            src/Indexer.cpp
            src/InverseSpaceTransform.cpp
            src/Lattice.cpp
            src/LatticeAssembler.cpp
            src/pointAutocorrelation.cpp
            src/SamplePointsGenerator.cpp 
            src/SparsePeakFinder.cpp
            )
 

add_library(indexer SHARED ${SOURCES})

set_target_properties(indexer PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS ON
)

find_package(Eigen3 3.3 NO_MODULE)
if(EIGEN3_FOUND)
    target_link_libraries (indexer Eigen3::Eigen)
else()
    message("Eigen not found, using internal Eigen.")
    include_directories(include/Eigen)
    
    install(
        DIRECTORY include/Eigen
        DESTINATION include
    )
endif()

install(TARGETS indexer 
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION lib
)

install(
    DIRECTORY include
    DESTINATION include
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "Eigen/*" EXCLUDE
    PATTERN "matlabLikeFunctions.h" EXCLUDE
    PATTERN "sortingByOtherValues.h" EXCLUDE
)

