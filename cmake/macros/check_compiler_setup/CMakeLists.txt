cmake_minimum_required(VERSION 3.13.4)
project(CheckCompilerSetup)
add_executable(CheckCompilerSetupExec dummy.cpp)

# Make sure some CUDA warning flags don't get deduplicated
string (REPLACE ";" " " TEST_COMPILE_FLAGS "${TEST_COMPILE_FLAGS}")
string(REGEX REPLACE "(-Xcudafe --diag_suppress=[^ ]+)" "\"SHELL:\\1\"" TEST_COMPILE_FLAGS "${TEST_COMPILE_FLAGS}")
separate_arguments(TEST_COMPILE_FLAGS UNIX_COMMAND "${TEST_COMPILE_FLAGS}")

target_compile_options(CheckCompilerSetupExec PRIVATE ${TEST_COMPILE_FLAGS})
set_property(TARGET CheckCompilerSetupExec PROPERTY LINK_FLAGS "${TEST_LINK_OPTIONS}")
target_link_libraries(CheckCompilerSetupExec ${TEST_LINK_LIBRARIES})

