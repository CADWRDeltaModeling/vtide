cmake_minimum_required(VERSION 3.0.2)
SET (CMAKE_Fortran_COMPILER_FORCED true)
SET (CMAKE_C_COMPILER_FORCED true)
project (VTide Fortran C)

include(FindLAPACK)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -extend-source 132")
endif()

add_library(versatile_tide analyze.for astro_time.for constituents.for file_io.for predict.for svd.for write_output.for)

add_executable(vt_analyze analyze_driver.for)
target_link_libraries(vt_analyze versatile_tide ${LAPACK_LIBRARIES})

add_executable(vt_predict predict_driver.for)
target_link_libraries(vt_predict versatile_tide ${LAPACK_LIBRARIES})