file(REMOVE_RECURSE
  "cuda"
  "cuda.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/cuda.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
