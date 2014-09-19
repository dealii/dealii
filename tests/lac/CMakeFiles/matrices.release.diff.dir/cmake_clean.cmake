FILE(REMOVE_RECURSE
  "CMakeFiles/matrices.release.diff"
  "matrices.release/diff"
  "matrices.release/output"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/matrices.release.diff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
