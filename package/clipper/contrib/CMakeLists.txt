project(clipper-contrib)

add_library(clipper-contrib STATIC 
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/convolution_search.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/edcalc.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/fffear.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/function_object_bases.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/mapfilter.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/originmatch.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/sfcalc.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/sfcalc_obs.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/sfscale.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/sfweight.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/skeleton.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/contrib/test_contrib.cpp
)

#target_include_directories(clipper-contrib PRIVATE ../.. ${WRK_DIR}/checkout/clipper ${WRK_DIR}/fftw ${WRK_DIR}/rfftw)

target_include_directories(clipper-contrib PRIVATE
        ${clipper-dependencies_SOURCE_DIR}/clipper
        ${fftw-dependencies_SOURCE_DIR}
)

set(clipper-contrib_HEADERS
${WRK_DIR}/checkout/clipper//clipper/contrib/originmatch.h
${WRK_DIR}/checkout/clipper//clipper/contrib/fffear.h
${WRK_DIR}/checkout/clipper//clipper/contrib/function_object_bases.h
${WRK_DIR}/checkout/clipper//clipper/contrib/test_contrib.h
${WRK_DIR}/checkout/clipper//clipper/contrib/sfscale.h
${WRK_DIR}/checkout/clipper//clipper/contrib/skeleton.h
${WRK_DIR}/checkout/clipper//clipper/contrib/sfcalc_obs.h
${WRK_DIR}/checkout/clipper//clipper/contrib/edcalc.h
${WRK_DIR}/checkout/clipper//clipper/contrib/sfweight.h
${WRK_DIR}/checkout/clipper//clipper/contrib/convolution_search.h
${WRK_DIR}/checkout/clipper//clipper/contrib/mapfilter.h
${WRK_DIR}/checkout/clipper//clipper/contrib/sfcalc.h
${WRK_DIR}/checkout/clipper//clipper/clipper-contrib.h
)

target_compile_options(clipper-contrib PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-contrib PROPERTIES PUBLIC_HEADER "${clipper-contrib_HEADERS}")

install(TARGETS clipper-contrib
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/contrib
)
