FILE(REMOVE_RECURSE
  "CMakeFiles/solver.debug.diff"
  "solver.debug/diff"
  "solver.debug/output"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/solver.debug.diff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
