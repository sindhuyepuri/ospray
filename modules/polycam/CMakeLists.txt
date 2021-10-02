


cmake_dependent_option(
    OSPRAY_MODULE_POLYCAM
    "Build a polynomial lens model"
    OFF
    OSPRAY_ENABLE_MODULES
    OFF
)

if (NOT OSPRAY_MODULE_POLYCAM)
    return()
endif()

add_subdirectory(ospray)