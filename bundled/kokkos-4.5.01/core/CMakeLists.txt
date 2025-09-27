if(NOT Kokkos_INSTALL_TESTING)
  add_subdirectory(src)
endif()

function(kokkos_add_benchmark_directory DIR_NAME)
  if(NOT Kokkos_ENABLE_BENCHMARKS)
    return()
  endif()

  add_subdirectory(${DIR_NAME})
endfunction()

kokkos_add_test_directories(unit_test)
kokkos_add_benchmark_directory(perf_test)
