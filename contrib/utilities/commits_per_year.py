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

"""Determine the number of non-merge commits to deal.II in each full year.

Execute this program from any directory within a git tree, including from the
top-level directory.

Usage:
        python contrib/utilities/commits_per_year.py
"""

import datetime
import subprocess

import matplotlib.pyplot as plt

commits_per_year = []

current_year = datetime.datetime.now().year
years = range(1997, current_year + 1)
for year in years:
    commits = subprocess.run(
        [
            "git",
            "log",
            "--no-merges",
            "--since",
            f"{year}/01/01",
            "--until",
            f"{year}/12/31",
            "--format=%H",
        ],
        stdout=subprocess.PIPE,
        check=True,
        text=True,
    ).stdout.splitlines()

    print(year, len(commits))
    commits_per_year.append(len(commits))

# Show the last bar in gray to indicate that the year isn't over yet:
bar_colors = ["C0"] * len(years)
bar_colors[-1] = "gray"

fig, ax = plt.subplots()
ax.bar(years, commits_per_year, color=bar_colors)
ax.set_xlabel("Year")
ax.set_ylabel("Number of non-merge commits")
ax.set_title("deal.II non-merge commits in a given year")
ax.set_xticks(years[::5])
fig.tight_layout()
plt.show()
