
IF(DEAL_II_INSTALL_EXAMPLES)
  #
  # Add a rule on how to install the examples:
  #

  SET(deal_ii_examples
    step-1 step-2 step-3 step-4 step-5 step-6 step-7 step-8 step-9
    step-10 step-11 step-12 step-13 step-14 step-15 step-16 step-17
    step-18 step-19 step-20 step-21 step-22 step-23 step-24 step-25
    step-26 step-27 step-28 step-29 step-30 step-31 step-32 step-33
    step-34 step-35 step-36 step-37 step-38 step-39 step-40 step-41
    step-42 step-43 step-44 step-45 step-46 step-47 step-48
    )

  FOREACH(step ${deal_ii_examples})
    INSTALL(DIRECTORY ${step}
      DESTINATION ${DEAL_II_EXAMPLES_RELDIR}
      COMPONENT examples
      PATTERN "doc" EXCLUDE
      )
  ENDFOREACH()

ENDIF()
