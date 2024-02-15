cmake_minimum_required(VERSION 3.12)
project(SAILS LANGUAGES C CXX)

SET(SAILS_CC_FLAGS "-g -O0 -w")
SET(SAILS_CXX_FLAGS "-g -O0 -w")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SAILS_CC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SAILS_CXX_FLAGS} ")

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

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
#
#add_subdirectory(clipper/gemmi)


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

add_library(
        sails_lib
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-find.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-lib.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-model.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-refine.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-score.cpp
)


add_executable(csails ${CMAKE_CURRENT_SOURCE_DIR}/src/sails.cpp)

set_target_properties(csails PROPERTIES OUTPUT_NAME sails)

target_include_directories(csails PUBLIC checkout/gemmi)

target_link_directories(csails
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

target_link_libraries(
        csails
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
        #        clipper-gemmi
        clipper-ccp4
)

set(SAILS_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-find.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-lib.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-model.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-refine.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sails-score.h
)

set_target_properties(csails PROPERTIES PUBLIC_HEADER "${SAILS_HEADERS}")

install(TARGETS csails
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/sails
)

# Detect the installed nanobind package and import it into CMake
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)

message (STATUS "NBDIR = ${NB_DIR}")
message (STATUS "${Python_EXECUTABLE} -m nanobind --cmake_dir")

list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)
nanobind_add_module(pysails ${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/sails-bind.cpp)

