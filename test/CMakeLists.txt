
message(STATUS "COMPILER_ID: ${CMAKE_CXX_COMPILER_ID} and SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
if (CMAKE_SYSTEM_NAME STREQUAL "Windows" AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  message(STATUS "Intel Googletest workaround applied")  
  add_compile_options(/w)
endif()
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(hs071 src/hs071_main.cpp src/hs071_nlp.cpp)
target_link_libraries(hs071 PRIVATE ${COINIPOPTLIBS} GTest::gtest_main)
target_include_directories(hs071 PRIVATE include
                                        ${IPOPT_PREFIX}/include/coin-or)
target_link_directories(hs071 PRIVATE ${IPOPT_PREFIX}/lib)
if (WIN32 AND (MINGW OR CMAKE_C_COMPILER_ID STREQUAL "GNU"))
    target_link_options(hs071 PRIVATE
        -static
        -static-libgcc
        -static-libstdc++)
endif()


include(GoogleTest)
gtest_discover_tests(hs071)