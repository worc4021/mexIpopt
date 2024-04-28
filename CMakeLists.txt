set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake CACHE STRING "Vcpkg toolchain file")
cmake_minimum_required(VERSION 3.20)
project(ipopt_mex LANGUAGES CXX)

find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM)

find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)

set(MKL_LINK "static")
set(MKL_THREADING "sequential")
set(MKL_INTERFACE_FULL "intel_lp64")
find_package(MKL REQUIRED)

set(COINIPOPTLIBS)
if (WIN32)
    set(CMAKE_CXX_FLAGS "/Wall /EHsc")
    list(APPEND COINIPOPTLIBS libipopt.lib libcoinhsl.lib libcoinmumps.lib libcoinmetis.lib)
else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra")
    list(APPEND COINIPOPTLIBS libipopt.a libcoinhsl.a libcoinmumps.a libcoinmetis.a)
endif()

add_subdirectory(MexUtilities)

matlab_add_mex(
    NAME ipopt_mex
    SRC src/main.cpp
    OUTPUT_NAME ipopt
    LINK_TO 
        MexUtilities::MexUtilities 
        ${COINIPOPTLIBS}
        MKL::MKL
    R2018a
)

matlab_add_mex(
    NAME model_eval_mex
    SRC src/modelEval.cpp
    OUTPUT_NAME model_eval
    LINK_TO 
        MexUtilities::MexUtilities 
        ${COINIPOPTLIBS}
        MKL::MKL
    R2018a
)

target_include_directories(ipopt_mex PRIVATE 
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${ipopt_INCLUDE_DIRS})

target_link_directories(ipopt_mex PRIVATE ${ipopt_LIBRARY_DIRS})

target_include_directories(model_eval_mex PRIVATE 
                            ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${ipopt_INCLUDE_DIRS})

target_link_directories(model_eval_mex PRIVATE ${ipopt_LIBRARY_DIRS})

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