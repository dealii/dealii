#!/usr/bin/env python3
## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2018 - 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

"""Plot the code-size history produced by count_lines.sh.

Run this script from the repository root after creating
``code-size-over-time.dat`` with ``count_lines.sh``:

    contrib/utilities/plot_code_size.py code-size-over-time.dat

The input contains whitespace-separated ``date source-lines test-lines``
records. The plot is displayed in an interactive window; it is not written to
a file. ``RELEASES`` defines the release arrows as version, release date,
arrow start and end heights, label height, and arrow line width.
"""

import argparse
from datetime import date
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# Each entry is version, date, arrow start/end, label height, and line width.
RELEASES = (
    ("3.0", "2000-04-06", 250000, 150000, 265000, 2),
    ("4.0", "2003-06-04", 360000, 230000, 375000, 2),
    ("5.0", "2004-05-21", 300000, 250000, 315000, 2),
    ("6.0", "2007-09-06", 430000, 330000, 445000, 2),
    ("7.0", "2011-01-09", 100000, 400000, 75000, 2),
    ("8.0", "2013-07-24", 150000, 450000, 125000, 2),
    ("9.0", "2018-05-11", 250000, 550000, 225000, 2),
    ("9.1", "2019-05-21", 300000, 600000, 275000, 1),
    ("9.2", "2020-05-20", 350000, 650000, 325000, 1),
    ("9.3", "2021-06-18", 400000, 700000, 375000, 1),
    ("9.4", "2022-06-24", 450000, 750000, 425000, 1),
    ("9.5", "2023-07-07", 500000, 800000, 475000, 1),
    ("9.6", "2024-08-11", 550000, 850000, 525000, 1),
    ("9.7", "2025-07-22", 600000, 900000, 575000, 1),
    ("9.8", "2026-08-09", 650000, 950000, 625000, 1),
)


def read_data(filename: Path) -> tuple[list[date], list[int], list[int]]:
    dates: list[date] = []
    source_lines: list[int] = []
    test_lines: list[int] = []

    with filename.open(encoding="utf-8") as data:
        for line_number, line in enumerate(data, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 3:
                raise ValueError(
                    f"{filename}:{line_number}: expected date, source lines, and test lines"
                )

            try:
                dates.append(date.fromisoformat(fields[0]))
                source_lines.append(int(fields[1]))
                test_lines.append(int(fields[2]))
            except ValueError as error:
                raise ValueError(f"{filename}:{line_number}: invalid data") from error

    if not dates:
        raise ValueError(f"{filename}: contains no data")

    # count_lines.sh traverses history backwards, but plots need chronological data.
    sorted_data = sorted(zip(dates, source_lines, test_lines))
    return (
        [entry[0] for entry in sorted_data],
        [entry[1] for entry in sorted_data],
        [entry[2] for entry in sorted_data],
    )


def main() -> None:
    # Keep the data file configurable while preserving the usual repository-root
    # invocation.
    parser = argparse.ArgumentParser(
        description="Plot the code-size history produced by count_lines.sh."
    )
    parser.add_argument(
        "input",
        type=Path,
        nargs="?",
        default=Path("code-size-over-time.dat"),
        help="input data file (default: %(default)s)",
    )
    arguments = parser.parse_args()

    dates, source_lines, test_lines = read_data(arguments.input)

    figure, axes = plt.subplots()
    axes.plot(
        dates,
        source_lines,
        color="tab:orange",
        label="Lines of code in source files",
    )
    axes.plot(dates, test_lines, color="tab:blue", label="Lines of code in tests")
    # A fixed reference rate makes long-term growth easier to compare.
    axes.plot(
        (date(2008, 1, 1), date(2026, 12, 31)),
        (290000, 1050000),
        color="gray",
        label="40,000 lines per year",
    )

    # The arrow head is at ``end``; early releases therefore point down while
    # later ones point up from below the curves.
    for version, release_date, start, end, label_position, line_width in RELEASES:
        release = date.fromisoformat(release_date)
        axes.annotate(
            "",
            xy=(release, end),
            xytext=(release, start),
            arrowprops={"arrowstyle": "->", "linewidth": line_width},
        )
        axes.annotate(version, xy=(release, label_position), ha="center")

    # Anchor five-year ticks at 1997 rather than at years divisible by five.
    axes.set_xlim(date(1997, 1, 1), date(2026, 12, 31))
    axes.set_xlabel("Date")
    axes.set_ylabel("Lines of code")
    axes.set_xticks([date(year, 1, 1) for year in range(1997, 2027, 5)])
    axes.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axes.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.0f}"))
    axes.legend(loc="upper left")
    plt.show()


if __name__ == "__main__":
    main()
