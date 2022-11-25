cmake_minimum_required(VERSION 3.3.0)
include(../setup_testsubproject.cmake)
project(testsuite CXX)

if(DEAL_II_WITH_TRILINOS AND DEAL_II_TRILINOS_WITH_BELOS)
  deal_ii_pickup_tests()
endif()
