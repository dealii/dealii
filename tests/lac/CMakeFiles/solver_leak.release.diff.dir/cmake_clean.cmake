FILE(REMOVE_RECURSE
  "CMakeFiles/solver_leak.release.diff"
  "solver_leak.release/diff"
  "solver_leak.release/output"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/solver_leak.release.diff.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
