FILE(REMOVE_RECURSE
  "CMakeFiles/solver_leak.debug.diff"
  "solver_leak.debug/diff"
  "solver_leak.debug/output"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/solver_leak.debug.diff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
