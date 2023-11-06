cmake_minimum_required(VERSION 3.2...3.5)
project(MiniCam LANGUAGES CXX)

# write .gitignore to CMAKE_BINARY_DIR
file(WRITE ${CMAKE_BINARY_DIR}/.gitignore "*")

add_library(minicam INTERFACE)
target_include_directories(minicam INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_executable(example example.cpp)
target_link_libraries(example minicam)
target_compile_options(example PRIVATE "$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

find_package(GTest QUIET)
if(GTest_FOUND)
    add_executable(unittest unittest.cpp)
    target_link_libraries(unittest PRIVATE GTest::gtest GTest::gtest_main)
endif()