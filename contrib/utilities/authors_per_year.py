#!/usr/bin/env python

## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

"""This script is used to determine the number of unique authors who
have contributed to deal.II in each (full) year since the beginning
of the project. Execute this program from any directory within a git tree,
including from the top level directory.

Usage:
        python contrib/utilities/authors_per_year.py
"""

import datetime
import subprocess
import matplotlib.pyplot as plt

authors_per_year = []

current_year = datetime.datetime.now().year
years = range(1997, current_year + 1)
for year in years:
    # Ask git for the logs from that particular year:
    log = subprocess.run(
        ["git", "log", "--since", f"{year}/01/01", "--until", f"{year}/12/31"],
        stdout=subprocess.PIPE,
    ).stdout.decode("cp437")

    # Then split off all lines that start with Author: and strip email
    # addresses:
    authors = []
    for line in log.splitlines():
        if "Author: " in line:
            author = line[8 : line.find(" <")]
            authors.append(author)

    # Make the list unique and show it:
    authors = list(set(authors))
    print(year, len(authors))

    authors_per_year.append(len(authors))

# Show the last bar in gray to indicate that the year isn't over yet:
bar_colors = ["C0"] * len(years)
bar_colors[-1] = "gray"

fig, ax = plt.subplots()
ax.bar(years, authors_per_year, color=bar_colors)
ax.set_xlabel("Year")
ax.set_ylabel("Number of unique authors")
ax.set_title("deal.II authors in a given year")
ax.set_xticks(years[::5])
fig.tight_layout()
plt.show()
