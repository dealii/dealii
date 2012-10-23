#
# A minimalistic CMakeLists.txt example for a small project using the
# deal.II library:
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 REQUIRED
  HINTS
    ../ ../../ ../../../ ../../../../../
    $ENV{DEAL_II_DIR} $ENV{DEAL_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(step-2)

ADD_EXECUTABLE(step-2 step-2.cc)
DEAL_II_SETUP_TARGET(step-2)

