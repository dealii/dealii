set(TARGET "step-40")

set(TARGET_SRC ${TARGET}.cc)

# Usually, you will not need to modify anything beyond this point...

cmake_minimum_required(VERSION 3.13.4)

find_package(deal.II 9.6.0 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

deal_ii_initialize_cached_variables()
project(${TARGET})
deal_ii_invoke_autopilot()
