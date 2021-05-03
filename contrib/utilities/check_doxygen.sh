#!/bin/bash

if [ ! -f contrib/utilities/checkdoxygen.py ]; then
  echo "*** This script must be run from the top-level directory of deal.II."
  exit 1
fi


find doc examples include \( -name "*.h" -o -name "*.dox" \) -print0 | xargs -0 -n 1 contrib/utilities/checkdoxygen.py
