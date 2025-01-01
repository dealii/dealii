# use a "Custom" build type to avoid auto populated compiler flags
set(CMAKE_BUILD_TYPE "Custom")
set(CMAKE_CONFIGURATION_TYPES "DEBUG;RELEASE")
cmake_minimum_required(VERSION 3.13.4)
project(CheckCompilerSetup)
add_executable(CheckCompilerSetupExec dummy.cpp)

target_compile_options(CheckCompilerSetupExec PRIVATE ${TEST_COMPILE_OPTIONS})
target_link_options(CheckCompilerSetupExec PRIVATE ${TEST_LINK_OPTIONS})
target_link_libraries(CheckCompilerSetupExec ${TEST_LINK_LIBRARIES})
