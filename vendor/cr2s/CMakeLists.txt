cmake_minimum_required(VERSION 3.1)
project(cR2S)

set(cR2S_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})


# Build
# =====

set(cR2S_HEADERS
    ${cR2S_INCLUDE_DIR}/cR2smodules.h
)

add_library(cR2S INTERFACE)
target_include_directories(Cr2S INTERFACE $<BUILD_INTERFACE:${cR2s_INCLUDE_DIR}>
                                         $<INSTALL_INTERFACE:include>)
# Installation
# ============

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS cR2S
        EXPORT ${PROJECT_NAME}-targets)

# Makes the project importable from the build directory
export(EXPORT ${PROJECT_NAME}-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")

install(FILES ${cR2S_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cR2S)

set(cR2S_CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for cR2SConfig.cmake")

configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION ${cR2S_CMAKECONFIG_INSTALL_DIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${cR2S_CMAKECONFIG_INSTALL_DIR})
install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${cR2S_CMAKECONFIG_INSTALL_DIR})
