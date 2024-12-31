set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake CACHE STRING "Vcpkg toolchain file")
cmake_minimum_required(VERSION 3.20)
project(ipopt_mex LANGUAGES CXX)

find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM)

option(IPOPT_PREFIX "Path to the IPOPT installation" "")
if (NOT IPOPT_PREFIX)
    message(FATAL_ERROR "IPOPT_PREFIX must be set to the path of the IPOPT installation")
else()
    if (WIN32)
        set(ENV{PKG_CONFIG_PATH} "${IPOPT_PREFIX}/lib/pkgconfig;$ENV{PKG_CONFIG_PATH}")
    else()
        set(ENV{PKG_CONFIG_PATH} "${IPOPT_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)

set(COINIPOPTLIBS)
if (WIN32)
    set(MKL_LINK "static")
    set(MKL_THREADING "sequential")
    set(MKL_INTERFACE_FULL "intel_lp64")
    find_package(MKL REQUIRED)
    set(CMAKE_CXX_FLAGS "/Wall /EHsc")
    list(APPEND COINIPOPTLIBS ${ipopt_LIBRARY_DIRS}/ipopt.lib ${ipopt_LIBRARY_DIRS}/coinhsl.lib ${ipopt_LIBRARY_DIRS}/coinmumps.lib ${ipopt_LIBRARY_DIRS}/libcoinmetis.lib MKL::MKL)
else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra")
    list(APPEND COINIPOPTLIBS PkgConfig::ipopt)
endif()

add_subdirectory(MexUtilities)

matlab_add_mex(
    NAME ipopt_mex
    SRC src/main.cpp
    OUTPUT_NAME ipopt
    LINK_TO 
        MexUtilities::MexUtilities 
        ${COINIPOPTLIBS}
    R2018a
)

matlab_add_mex(
    NAME model_eval_mex
    SRC src/modelEval.cpp
    OUTPUT_NAME model_eval
    LINK_TO 
        MexUtilities::MexUtilities 
        ${COINIPOPTLIBS}
    R2018a
)

target_include_directories(ipopt_mex PRIVATE 
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${ipopt_INCLUDE_DIRS})

set_target_properties(ipopt_mex PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

target_include_directories(model_eval_mex PRIVATE 
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${ipopt_INCLUDE_DIRS})
set_target_properties(model_eval_mex PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")                          

message(STATUS "binary path: ${CMAKE_CURRENT_BINARY_DIR}")

include(CTest)
matlab_add_unit_test(NAME HS71_DirectInterface
                     UNITTEST_FILE hs71Direct.m
                     ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/matlab $<TARGET_FILE_DIR:ipopt_mex>)

matlab_add_unit_test(NAME HS71_ModelInterface
                     UNITTEST_FILE hs71Model.m`
                     ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/matlab $<TARGET_FILE_DIR:ipopt_mex>)

matlab_add_unit_test(NAME Himmelblau_DirectInterface
                     UNITTEST_FILE HimmelblauDirect.m
                     ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/matlab $<TARGET_FILE_DIR:ipopt_mex>)

matlab_add_unit_test(NAME Himmelblau_ModelInterface
                     UNITTEST_FILE HimmelblauModel.m
                     ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/matlab $<TARGET_FILE_DIR:ipopt_mex>)