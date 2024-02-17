cmake_minimum_required(VERSION 3.15...3.26)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

project(sails LANGUAGES CXX)

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
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()

add_subdirectory(fftw)
add_subdirectory(rfftw)
add_subdirectory(ccp4)
add_subdirectory(mmdb2)
add_subdirectory(gemmi)
add_subdirectory(clipper/minimol)
add_subdirectory(clipper/core)
add_subdirectory(clipper/contrib)
add_subdirectory(clipper/cns)
add_subdirectory(clipper/mmdb)
add_subdirectory(clipper/cif)
add_subdirectory(clipper/phs)
add_subdirectory(clipper/ccp4)
add_subdirectory(clipper/gemmi)

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        checkout/mmdb2
        checkout/clipper
        checkout/libccp4/ccp4
        checkout/gemmi/include
)

link_directories(
        mmdb2
        clipper
        ccp4
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Try to import all Python components potentially needed by nanobind
find_package(Python 3.8
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
        src/python/py-sails-input.cpp
)

target_include_directories(sails_module PUBLIC
        checkout/gemmi
        mmdb2
        clipper/mmdb
        clipper/cif
        clipper/core
        clipper/contrib
        clipper/ccp4
        clipper/gemmi
        clipper/minimol
        fftw
        rfftw
        ccp4
)


target_link_directories(sails_module
        PUBLIC
        mmdb2
        clipper/mmdb
        clipper/cif
        clipper/core
        clipper/contrib
        clipper/ccp4
        clipper/gemmi
        clipper/minimol
        fftw
        rfftw
        ccp4
)

add_library(
        sails_lib
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-find.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-lib.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-model.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-refine.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-score.cpp
)

target_link_libraries(
        sails_module
        PRIVATE
        gemmi_cpp
        sails_lib
        ccp4
        fftw
        rfftw
        mmdb2
        clipper-mmdb
        clipper-cif
        clipper-core
        clipper-contrib
        clipper-minimol
        clipper-ccp4
)
# Install directive for scikit-build-core
install(TARGETS sails_module LIBRARY DESTINATION pysails)

#add_executable(csails ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/csails.cpp)
#set_target_properties(csails PROPERTIES OUTPUT_NAME sails)
#target_include_directories(csails PUBLIC checkout/gemmi)
#target_link_directories(csails
#        PUBLIC
#        mmdb2
#        clipper/mmdb
#        clipper/cif
#        clipper/core
#        clipper/contrib
#        clipper/ccp4
#        clipper/gemmi
#        clipper/minimol
#        fftw
#        rfftw
#        ccp4
#)
#target_link_libraries(
#        csails
#        PRIVATE
#        gemmi_cpp
#        sails_lib
#        ccp4
#        fftw
#        rfftw
#        mmdb2
#        clipper-mmdb
#        clipper-cif
#        clipper-core
#        clipper-contrib
#        clipper-minimol
#        clipper-ccp4
#)
#set(SAILS_HEADERS
#        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-find.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-lib.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-model.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-refine.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/sails-score.h
#)
#set_target_properties(csails PROPERTIES PUBLIC_HEADER "${SAILS_HEADERS}")
#install(TARGETS csails
#        LIBRARY DESTINATION lib
#        PUBLIC_HEADER DESTINATION include/sails)
