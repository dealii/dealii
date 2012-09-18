
MACRO(STRIP_FLAG variable flag)
  STRING(REGEX REPLACE " ${flag}" "" ${variable} ${${variable}})
ENDMACRO()

