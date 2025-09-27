# Example Testing

This folder, that is part of the deal.II testsuite, contains the logic
to run the deal.II examples (located under ./examples/ in the
repository) as automated tests using CTest. Instead of storing copies
of the source code of these example, we only store a fairly small
textual diff for each example that, for example, limit the number of
time steps or number of mesh refinement steps. Wanting to test the
current form of the tutorial programs, perhaps up to these small
modifications, is motivated by the lesson that
copies of examples very quickly become out of date and we repeatedly
fail to detect if an example program breaks due to a change in the
library.

## Usage

When running the ``setup_tests`` or ``setup_tests_examples`` target in
a deal.II build directory, the diff files will be used to create
temporary .cc files located at
``build/tests/examples/source/step-*.cc``.

For development of new examples, we suggest first creating an empty
.diff and .output file in the tests/examples/ folders and then running
``setup_tests`` to pick up the new test.

To modify the source code, go into ``build/tests/examples/source/``
and edit the .cc file there. Afterwards, you can update the .diff
files by running the ``update_diffs`` target in the folder
``build/tests/examples/``, which runs the ``update_diffs.sh`` script
to update the diffs.
