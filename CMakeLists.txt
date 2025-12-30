cmake_minimum_required(VERSION 3.25)
project(ipopt_mex LANGUAGES CXX)

find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")


if (DEFINED ENV{BUILD_BUILDNUMBER})
    set(BUILD_NUMBER $ENV{BUILD_BUILDNUMBER})
else()
    set(BUILD_NUMBER "local")
endif()

find_package(Git)
if(Git_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE IPOPT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "Git commit hash: ${IPOPT_HASH}")
else()
    set(IPOPT_HASH "unknown version")
    message(WARNING "Git not found, cannot determine commit hash")
endif()

set(COINIPOPTLIBS "")
set(IPOPT_PREFIX "C:/repos/Ipopt/build" CACHE PATH "Path to the IPOPT installation")
if (IS_DIRECTORY "${IPOPT_PREFIX}")
    if (WIN32)
        set(MKL_LINK "static")
        set(MKL_THREADING "sequential")
        set(MKL_INTERFACE_FULL "intel_lp64")
        find_package(MKL REQUIRED)
        set(CMAKE_CXX_FLAGS "/Wall /EHsc")
        list(APPEND COINIPOPTLIBS   ${IPOPT_PREFIX}/lib/ipopt.lib 
                                    ${IPOPT_PREFIX}/lib/coinhsl.lib 
                                    ${IPOPT_PREFIX}/lib/coinmumps.lib 
                                    ${IPOPT_PREFIX}/lib/libcoinmetis.lib
                                    MKL::MKL)
    else(WIN32)
        set(CMAKE_CXX_FLAGS "-Wall -Wextra -fpie")
        set(ENV{PKG_CONFIG_PATH} "${IPOPT_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)
        list(APPEND COINIPOPTLIBS PkgConfig::ipopt)
    endif(WIN32)
else()
    message(FATAL_ERROR "IPOPT_PREFIX must be set to the path of the IPOPT installation, found '${IPOPT_PREFIX}'")
endif()
configure_file(version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/version.h @ONLY)

include(FetchContent)
set(USE_EIGEN OFF)
set(BUILD_TESTS OFF)
FetchContent_Declare(
    MexUtilities
    GIT_REPOSITORY https://github.com/worc4021/MexUtilities
    GIT_TAG        main
    EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(MexUtilities)

matlab_add_mex(
    NAME ipopt_mex
    SRC src/main.cpp
    OUTPUT_NAME ipopt
    LINK_TO 
        MexUtilities::MexUtilities 
        ${COINIPOPTLIBS}
    R2018a
)

target_include_directories(ipopt_mex PRIVATE 
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${CMAKE_CURRENT_BINARY_DIR}/include
                            ${IPOPT_PREFIX}/include/coin-or)

target_link_directories(ipopt_mex PRIVATE ${IPOPT_PREFIX}/lib)


message(STATUS "binary path: ${CMAKE_CURRENT_BINARY_DIR}")

enable_testing()
include(CTest)
add_subdirectory(test)

matlab_add_unit_test(NAME HS71_DirectInterface
                     UNITTEST_FILE hs71Direct.m
                     ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/matlab $<TARGET_FILE_DIR:ipopt_mex>)

matlab_add_unit_test(NAME Himmelblau_DirectInterface
                     UNITTEST_FILE HimmelblauDirect.m
                     ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/matlab $<TARGET_FILE_DIR:ipopt_mex>)

matlab_add_unit_test(NAME OptionsTest
                    UNITTEST_FILE OptionsTest.m
                    ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/matlab $<TARGET_FILE_DIR:ipopt_mex>)