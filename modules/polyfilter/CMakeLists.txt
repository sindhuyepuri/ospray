
cmake_dependent_option(
    OSPRAY_MODULE_POLYFILTER
    "Build filter that will simulate abberrations attributed to achromatic doublet lens"
    OFF
    OSPRAY_ENABLE_MODULES
    OFF
)

if (NOT OSPRAY_MODULE_POLYFILTER)
    return()
endif()

include(GenerateExportHeader)

add_library(ospray_module_polyfilter SHARED ${OSPRAY_RESOURCE} PolyFilterFrameOp.cpp)

generate_export_header(ospray_module_polyfilter)

target_link_libraries(ospray_module_polyfilter
PRIVATE
    ospray_module_ispc
)

target_include_directories(ospray_module_polyfilter
PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

ospray_install_library(ospray_module_polyfilter lib)