FILE(REMOVE_RECURSE
  "CMakeFiles/trace.debug.diff"
  "trace.debug/diff"
  "trace.debug/output"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/trace.debug.diff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
