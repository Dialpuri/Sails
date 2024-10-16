cmake_minimum_required(VERSION 3.15...3.26)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)
project(sails LANGUAGES CXX)
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (NOT SKBUILD)
    message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve -Ceditable.rebuild=true .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif ()

find_package(ZLIB REQUIRED)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(WRK_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

include(FetchContent)
# Download Clipper
FetchContent_Declare(
        clipper-dependencies
        URL http://www.ysbl.york.ac.uk/jsd523/clipper-2.1.20201109.tar.gz
)

# Download MMDB
FetchContent_Declare(
        mmdb2-dependencies
        URL http://www.ysbl.york.ac.uk/jsd523/mmdb2-2.0.22.tar.gz
)

# Download FFTW
FetchContent_Declare(
        fftw-dependencies
        URL http://www.fftw.org/fftw-2.1.5.tar.gz
)

# Download CCP4
FetchContent_Declare(
        ccp4-dependencies
        URL http://www.ysbl.york.ac.uk/jsd523/libccp4-8.0.0.tar.gz
)

# Download Gemmi
FetchContent_Declare(
        gemmi-dependencies
        URL http://www.ysbl.york.ac.uk/jsd523/gemmi-0.6.5.tar.gz
)

FetchContent_MakeAvailable(clipper-dependencies mmdb2-dependencies fftw-dependencies ccp4-dependencies gemmi-dependencies)

add_subdirectory(${WRK_DIR}/clipper/minimol)
add_subdirectory(${WRK_DIR}/clipper/core)
add_subdirectory(${WRK_DIR}/clipper/contrib)
add_subdirectory(${WRK_DIR}/clipper/cns)
add_subdirectory(${WRK_DIR}/clipper/mmdb)
add_subdirectory(${WRK_DIR}/clipper/cif)
add_subdirectory(${WRK_DIR}/clipper/phs)
add_subdirectory(${WRK_DIR}/clipper/ccp4)
add_subdirectory(${WRK_DIR}/fftw)
add_subdirectory(${WRK_DIR}/rfftw)
add_subdirectory(${WRK_DIR}/ccp4)
add_subdirectory(${WRK_DIR}/mmdb2)
add_subdirectory(${WRK_DIR}/gemmi)

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${clipper-dependencies_SOURCE_DIR}
        ${gemmi-dependencies_SOURCE_DIR}/include
        ${mmdb2-dependencies_SOURCE_DIR}
)

link_directories(
        ${CMAKE_BINARY_DIR}/mmdb2
        ${CMAKE_BINARY_DIR}/clipper
        ${CMAKE_BINARY_DIR}/ccp4
        ${WRK_DIR}/src
        ${CMAKE_BINARY_DIR}
)

add_library(
        sails_lib
        STATIC
        # Source
        ${WRK_DIR}/src/cpp/sails-json.cpp
        ${WRK_DIR}/src/cpp/sails-vector.cpp
        ${WRK_DIR}/src/cpp/sails-sequence.cpp
        ${WRK_DIR}/src/cpp/sails-topology.cpp
        ${WRK_DIR}/src/cpp/sails-glycan.cpp
        ${WRK_DIR}/src/cpp/sails-linkage.cpp
        ${WRK_DIR}/src/cpp/sails-utils.cpp
        ${WRK_DIR}/src/cpp/sails-refine.cpp
        ${WRK_DIR}/src/cpp/sails-cif.cpp
        ${WRK_DIR}/src/cpp/sails-dot.cpp
        ${WRK_DIR}/src/cpp/sails-telemetry.cpp
        ${WRK_DIR}/src/cpp/sails-solvent.cpp
        ${WRK_DIR}/src/cpp/sails-wurcs.cpp

        ${WRK_DIR}/src/cpp/sails-morph.cpp
        # Density
        ${WRK_DIR}/src/cpp/density/sails-density.cpp
        ${WRK_DIR}/src/cpp/density/sails-xtal-density.cpp
        ${WRK_DIR}/src/cpp/density/sails-em-density.cpp

        # SNFG
        ${WRK_DIR}/src/cpp/snfg/sails-snfg.cpp
        ${WRK_DIR}/src/cpp/snfg/sails-snfg-shape.cpp
        ${WRK_DIR}/src/cpp/snfg/shapes/circle.cpp
        ${WRK_DIR}/src/cpp/snfg/shapes/square.cpp
        ${WRK_DIR}/src/cpp/snfg/shapes/rectangle.cpp
        ${WRK_DIR}/src/cpp/snfg/shapes/triangle.cpp
        ${WRK_DIR}/src/cpp/snfg/shapes/diamond.cpp

        # Bindings
        ${WRK_DIR}/src/bindings/sails-gemmi-bindings.cpp

        # Third Party
        ${WRK_DIR}/src/third-party/simdjson.cpp
)

# target_include_directories(csails PUBLIC checkout/gemmi)
set(LINK_LIBS "gemmi_cpp sails_lib ccp4 fftw rfftw mmdb2 clipper-mmdb clipper-cif clipper-core clipper-contrib clipper-minimol clipper-ccp4 clipper-gemmi")

target_link_directories(sails_lib
        PUBLIC
        ${CMAKE_BINARY_DIR}/clipper/mmdb
        ${CMAKE_BINARY_DIR}/clipper/cif
        ${CMAKE_BINARY_DIR}/clipper/ccp4
        ${CMAKE_BINARY_DIR}/clipper/minimol
        ${CMAKE_BINARY_DIR}/clipper/core
        ${CMAKE_BINARY_DIR}/clipper/contrib
        ${CMAKE_BINARY_DIR}/fftw
        ${CMAKE_BINARY_DIR}/rfftw
        ${CMAKE_BINARY_DIR}/ccp4
        ${CMAKE_BINARY_DIR}/mmdb2
        ${CMAKE_BINARY_DIR}/gemmi
)

target_link_libraries(
        sails_lib
        PUBLIC
        clipper-mmdb
        clipper-cif
        clipper-ccp4
        clipper-minimol
        clipper-core
        clipper-contrib
        rfftw
        fftw
        ccp4
        mmdb2
        gemmi_cpp
)

file(GLOB sails_HEADERS ${WRK_DIR}/src/include/*.h)

if (BUILD_EXECUTABLE STREQUAL "ON")
    add_executable(csails ${WRK_DIR}/src/cpp/sails.cpp)
    set_target_properties(csails PROPERTIES OUTPUT_NAME sails_exec)
    target_compile_features(csails INTERFACE cxx_std_17)
    target_compile_options(csails PRIVATE -g)

    target_link_directories(csails
            PUBLIC
            ${CMAKE_BINARY_DIR}/clipper/mmdb
            ${CMAKE_BINARY_DIR}/clipper/cif
            ${CMAKE_BINARY_DIR}/clipper/ccp4
            ${CMAKE_BINARY_DIR}/clipper/minimol
            ${CMAKE_BINARY_DIR}/clipper/core
            ${CMAKE_BINARY_DIR}/clipper/contrib
            ${CMAKE_BINARY_DIR}/fftw
            ${CMAKE_BINARY_DIR}/rfftw
            ${CMAKE_BINARY_DIR}/ccp4
            ${CMAKE_BINARY_DIR}/mmdb2
            ${CMAKE_BINARY_DIR}/gemmi
    )

    target_link_libraries(
            csails
            PUBLIC
            sails_lib
            clipper-contrib
            clipper-mmdb
            clipper-cif
            clipper-ccp4
            clipper-minimol
            clipper-core
            rfftw
            fftw
            ccp4
            mmdb2
            gemmi_cpp
            ZLIB::ZLIB
    )

    set_target_properties(csails PROPERTIES PUBLIC_HEADER "${sails_HEADERS}")
    install(TARGETS csails
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION ${WRK_DIR}/include/sails)

endif ()

if (TESTING)
    include(FetchContent)
    FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    add_executable(
            sails_tests
            ${WRK_DIR}/src/cpp/tests/sails-glycan-test.cpp
    #        ${WRK_DIR}/src/cpp/tests/sails-linkage-test.cpp
    )
    target_link_libraries(
            sails_tests
            sails_lib
            ZLIB::ZLIB  #for writing tests
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(sails_tests)
endif()


#if (SKBUILD)
    # Try to import all Python components potentially needed by nanobind
    find_package(Python 3.10
            REQUIRED COMPONENTS Interpreter Development.Module
            OPTIONAL_COMPONENTS Development.SABIModule)

    # Import nanobind through CMake's find_package mechanism
    find_package(nanobind CONFIG REQUIRED)

    # We are now ready to compile the actual extension module
    nanobind_add_module(
            sails_module

            STABLE_ABI
            NB_STATIC

            # Python bound functions here
            ${WRK_DIR}/src/bindings/python_sails.cpp
    )


    target_include_directories(sails_module PUBLIC
            ${WRK_DIR}/checkout/gemmi
            ${WRK_DIR}/mmdb2
            ${WRK_DIR}/clipper/mmdb
            ${WRK_DIR}/clipper/cif
            ${WRK_DIR}/clipper/core
            ${WRK_DIR}/clipper/contrib
            ${WRK_DIR}/clipper/ccp4
            ${WRK_DIR}/clipper/minimol
            ${WRK_DIR}/fftw
            ${WRK_DIR}/rfftw
            ${WRK_DIR}/ccp4

    )


    target_link_directories(sails_module
            PUBLIC
            ${CMAKE_BINARY_DIR}/clipper/mmdb
            ${CMAKE_BINARY_DIR}/clipper/cif
            ${CMAKE_BINARY_DIR}/clipper/ccp4
            ${CMAKE_BINARY_DIR}/clipper/minimol
            ${CMAKE_BINARY_DIR}/clipper/contrib
            ${CMAKE_BINARY_DIR}/clipper/core
            ${CMAKE_BINARY_DIR}/fftw
            ${CMAKE_BINARY_DIR}/rfftw
            ${CMAKE_BINARY_DIR}/ccp4
            ${CMAKE_BINARY_DIR}/mmdb2
            ${CMAKE_BINARY_DIR}/gemmi
    )


    target_link_libraries(
            sails_module
            PRIVATE
            sails_lib
            clipper-mmdb
            clipper-cif
            clipper-ccp4
            clipper-minimol
            clipper-contrib
            clipper-core
            rfftw
            fftw
            ccp4
            mmdb2
            gemmi_cpp
            ZLIB::ZLIB

    )

    install(TARGETS sails_module LIBRARY DESTINATION sails)
#endif ()
