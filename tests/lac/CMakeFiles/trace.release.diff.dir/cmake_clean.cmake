FILE(REMOVE_RECURSE
  "CMakeFiles/trace.release.diff"
  "trace.release/diff"
  "trace.release/output"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/trace.release.diff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
