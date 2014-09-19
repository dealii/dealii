FILE(REMOVE_RECURSE
  "CMakeFiles/solver.release.diff"
  "solver.release/diff"
  "solver.release/output"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/solver.release.diff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
