cmake_minimum_required(VERSION 3.25)
project(ipopt_mex LANGUAGES CXX)


set(ARL_MATLAB_VERSION "R2022b" CACHE STRING "Matlab version to build against")

set(Matlab_ROOT_DIR "C:/Program Files/MATLAB/${ARL_MATLAB_VERSION}" CACHE PATH "Path to the MATLAB installation")
find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if (MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()


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
if (WIN32 AND (MINGW OR CMAKE_C_COMPILER_ID STREQUAL "GNU"))
    set(IPOPT_PREFIX_DEFAULT "C:/repos/Ipopt/build/mingw")
elseif (WIN32)
    set(IPOPT_PREFIX_DEFAULT "C:/repos/Ipopt/build/intel")
else()
    set(IPOPT_PREFIX_DEFAULT "C:/repos/Ipopt/build")
endif()
set(IPOPT_PREFIX "${IPOPT_PREFIX_DEFAULT}" CACHE PATH "Path to the IPOPT installation")
if (IS_DIRECTORY "${IPOPT_PREFIX}")
    if (WIN32 AND MSVC)
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
    elseif (WIN32 AND (MINGW OR CMAKE_C_COMPILER_ID STREQUAL "GNU"))
        set(CMAKE_CXX_FLAGS "-Wall -Wextra")
        set(ENV{PKG_CONFIG_PATH} "${IPOPT_PREFIX}/lib/pkgconfig")
        find_package(PkgConfig REQUIRED)
        # Force find_library to prefer static archives for the MinGW runtime
        # libraries pulled in by the ipopt pkg-config Requires.private chain
        # (libgfortran, libquadmath, libpthread, libwinpthread, ...).  Without
        # this, pkg_check_modules prefers the *.dll.a import libs and the mex
        # ends up depending on libgfortran-5.dll, libwinpthread-1.dll, etc.
        set(_orig_find_lib_suffixes "${CMAKE_FIND_LIBRARY_SUFFIXES}")
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)
        set(CMAKE_FIND_LIBRARY_SUFFIXES "${_orig_find_lib_suffixes}")
        list(APPEND COINIPOPTLIBS PkgConfig::ipopt)
    else()
        set(CMAKE_CXX_FLAGS "-Wall -Wextra -fpie")
        set(ENV{PKG_CONFIG_PATH} "${IPOPT_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)
        list(APPEND COINIPOPTLIBS PkgConfig::ipopt)
    endif()
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

if (WIN32 AND (MINGW OR CMAKE_C_COMPILER_ID STREQUAL "GNU"))
    target_link_options(ipopt_mex PRIVATE
        -static
        -static-libgcc
        -static-libstdc++)
endif()


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