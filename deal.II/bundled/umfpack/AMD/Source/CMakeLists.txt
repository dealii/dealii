SET(src_amd
  amd_1.c
  amd_2.c
  amd_aat.c
  amd_control.c
  amd_defaults.c
  amd_dump.c
  amd_info.c
  amd_order.c
  amd_postorder.c
  amd_post_tree.c
  amd_preprocess.c
  amd_valid.c
  )


DEAL_II_ADD_C_LIBRARY(obj_amd_int OBJECT ${src_amd})
DEAL_II_ADD_DEFINITIONS(obj_amd_int "DINT")

DEAL_II_ADD_C_LIBRARY(obj_amd_long OBJECT ${src_amd})
DEAL_II_ADD_DEFINITIONS(obj_amd_long "DLONG")

DEAL_II_ADD_C_LIBRARY(obj_amd_global OBJECT "amd_global.c")
