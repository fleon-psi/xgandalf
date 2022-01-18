cmake_minimum_required(VERSION 3.1)
project(xgandalf VERSION 1.0)

include(GNUInstallDirs)

option(XGANDALF_BUILD_EXECUTABLE "Build the test executable for xgandalf" OFF)
option(USE_INSTALLED_PRECOMUTED_DATA "Use the installation path for getting the precomputed data 
                                       (as oposite to the source location)" OFF)

 
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo")
endif()
 
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
 
# user defined options for MKL
option(MKL_USE_parallel "Use MKL parallel" False)
option(MKL_USE_sdl "Single Dynamic Library or static/dynamic" False)
set(BLA_STATIC True CACHE STRING "Use static linking instead of dynamic linking for MKL")
set(MKL_USE_interface "lp64" CACHE STRING "for Intel(R)64 compatible arch: ilp64/lp64 or for ia32 arch: cdecl/stdcall")


set(BOOST_ALL_NO_LIB 1)
find_package(Boost 1.60.0)
find_package(Eigen3 3.3.7 NO_MODULE)
find_package(MKL)

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
			src/refinement.cpp
			src/samplePointsFiltering.cpp
            src/SamplePointsGenerator.cpp 
            src/SparsePeakFinder.cpp

			src/ReciprocalToRealProjection.cpp
			src/SimpleMonochromaticDiffractionPatternPrediction.cpp
			src/SimpleMonochromaticProjection.cpp

			src/adaptions/crystfel/IndexerPlain.cpp
			src/adaptions/crystfel/ExperimentSettings.cpp
			src/adaptions/crystfel/indexerData.cpp
			src/adaptions/crystfel/lattice.cpp

			src/adaptions/crystfel/projectionData.cpp
			src/adaptions/crystfel/SimpleMonochromaticDiffractionPatternPrediction.cpp
  )

set(SOURCES_test src/main.cpp
				 src/tests.cpp
)

if(XGANDALF_BUILD_EXECUTABLE)
	add_executable(xgandalf ${SOURCES} ${SOURCES_test})
else(XGANDALF_BUILD_EXECUTABLE)
	add_library(xgandalf STATIC ${SOURCES})
endif(XGANDALF_BUILD_EXECUTABLE)

target_include_directories(xgandalf PUBLIC include)

if(USE_INSTALLED_PRECOMUTED_DATA)
	target_compile_definitions(xgandalf PRIVATE -DPRECOMPUTED_DATA_DIR=${CMAKE_INSTALL_FULL_DATADIR}/xgandalf/precomputedSamplePoints)
else(USE_INSTALLED_PRECOMUTED_DATA)
	target_compile_definitions(xgandalf PRIVATE -DPRECOMPUTED_DATA_DIR=${PROJECT_SOURCE_DIR}/precomputedSamplePoints)
endif(USE_INSTALLED_PRECOMUTED_DATA)

if(Boost_FOUND)
	target_include_directories(xgandalf SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
	message("Boost not found, using internal Boost.")
    target_include_directories(xgandalf SYSTEM PUBLIC include/Boost)

	install(
        DIRECTORY include/Boost/
        DESTINATION include/boost
    )
endif(Boost_FOUND)

if(EIGEN3_FOUND)
    target_link_libraries (xgandalf PUBLIC Eigen3::Eigen)
else()
    message("Eigen not found, using internal Eigen.")
    target_include_directories(xgandalf SYSTEM PUBLIC include/Eigen)
    
    install(
        DIRECTORY include/Eigen/
        DESTINATION include/eigen3
    )
endif()

if (MKL_FOUND)
	target_compile_definitions(xgandalf PUBLIC ${MKL_DEFINITIONS} USE_MKL EIGEN_USE_MKL_ALL)
	target_link_libraries(xgandalf PRIVATE ${MKL_LIBRARIES})
	target_include_directories(xgandalf PRIVATE ${MKL_INCLUDE_DIR})
endif(MKL_FOUND)

# Test whether the compiler is Microsoft Visual C(++).
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  target_compile_options(xgandalf PUBLIC /W2 /wd4305 /wd4244 /wd4099)
else(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  target_compile_options(xgandalf PRIVATE -Wall)
endif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
 

set_target_properties(xgandalf PROPERTIES 
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)



install(TARGETS xgandalf 
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
    DIRECTORY include/
    DESTINATION include/xgandalf
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "Eigen/*" EXCLUDE
	PATTERN "Eigen" EXCLUDE
	PATTERN "Boost" EXCLUDE
	PATTERN "Boost/*" EXCLUDE
)

install(
    DIRECTORY precomputedSamplePoints/
    DESTINATION share/xgandalf/precomputedSamplePoints
    FILES_MATCHING
    PATTERN "*"
	PATTERN "*.m" EXCLUDE
)

# xgandalf.pc
configure_file(xgandalf.pc.in xgandalf.pc)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/xgandalf.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)