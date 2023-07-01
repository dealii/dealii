cmake_minimum_required(VERSION 3.13.4)
include(../scripts/setup_testsubproject.cmake)
project(testsuite CXX)

# create tests from diffs:
execute_process(
    COMMAND "./update.sh"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

deal_ii_pickup_tests()
