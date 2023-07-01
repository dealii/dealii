#
# Glob all test summary files and remove them
#

if("${_summary_files_prefix}" STREQUAL "")
  message(FATAL_ERROR "The variable _summary_files_prefix is not set.")
endif()

file(GLOB _files "${_summary_files_prefix}_*")
if(NOT "${_files}" STREQUAL "")
  file(REMOVE ${_files})
endif()
