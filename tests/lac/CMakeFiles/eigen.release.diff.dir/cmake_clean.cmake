FILE(REMOVE_RECURSE
  "CMakeFiles/eigen.release.diff"
  "eigen.release/diff"
  "eigen.release/output"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/eigen.release.diff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
