## Copyright 2009-2019 Intel Corporation
## SPDX-License-Identifier: Apache-2.0

add_library(ospray_module_polycam SHARED)

ispc_include_directories(${RKCOMMON_INCLUDE_DIRS})

ispc_target_add_sources(ospray_module_polycam
  # the cpp file that contains all the plugin code - parsing
  # parameters in ospCommit(), creating and registering the object,
  # building accel structures, etc
  camera/PolynomialCamera.cpp

  # the ispc file that contains the plugins for all vector code - ie,
  # for ray-primitive intersection and 'postIntersect' (reporting info
  # on a previously computed ray-prim intersection)
  camera/PolynomialCamera.ispc

  # and finally, the module init code (not doing much, but must be there)

)

target_link_libraries(ospray_module_polycam PUBLIC ospray_module_ispc)
