#
# Glob all test summary files and print them
#

if("${_summary_files_prefix}" STREQUAL "")
  message(FATAL_ERROR "The variable _summary_files_prefix is not set.")
endif()

file(GLOB _files "${_summary_files_prefix}_*")
foreach(_file ${_files})
  file(STRINGS "${_file}" _summary_string REGEX "Test category ")
  if(NOT "${_summary_string}" STREQUAL "")
    message("${_summary_string}")
  endif()
endforeach()
