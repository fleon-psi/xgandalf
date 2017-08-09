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
 
# Test whether the compiler is Microsoft Visual C(++).
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W2 /wd4305 /wd4244 /wd4099")
endif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
 
set(BOOST_ALL_NO_LIB 1)
find_package(Boost 1.60.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(include)

set(SOURCES src/Dbscan.cpp
            src/DetectorToReciprocalSpaceTransform.cpp
            src/ExperimentSettings.cpp
            src/HillClimbingOptimizer.cpp
            src/IndexerAutocorrPrefit.cpp
			src/IndexerBase.cpp
			src/IndexerPlain.cpp
            src/InverseSpaceTransform.cpp
            src/Lattice.cpp
            src/LatticeAssembler.cpp
            src/pointAutocorrelation.cpp
			src/samplePointsFiltering.cpp
            src/SamplePointsGenerator.cpp 
            src/SparsePeakFinder.cpp
			src/tests.cpp

			src/adaptions/crystfel/IndexerPlain.cpp
			src/adaptions/crystfel/ExperimentSettings.cpp
			src/adaptions/crystfel/indexerData.cpp

			src/main.cpp
            )
 

add_library(indexer SHARED ${SOURCES})
#add_executable(indexer ${SOURCES})

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
    include_directories(SYSTEM include/Eigen)
    
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
)

