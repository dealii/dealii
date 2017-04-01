# Contributing to deal.II
deal.II is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive and
participatory community so we are happy that you are interested in
participating. This document partly duplicates some extra information on
participating available
[on the project website](https://www.dealii.org/participate.html).

## Getting started with git and GitHub
If you are new to using git or the GitHub platform you may find it helpful
to review
[lecture 32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html),
which should be enough to help you get started with using GitHub or
possibly contributing to the library itself.

## Asking and answering questions about the library
The deal.II community maintains an active mailing list hosted by google
groups [here](https://groups.google.com/forum/#!forum/dealii). The mailing
list is for questions about the library of all levels. For more information
see the website page describing our mailing lists
[here](https://www.dealii.org/mail.html).

## Bug reports
![bug report](https://cdn.meme.am/instances/500x/61103043/liam-neeson-taken-i-will-find-the-bug-and-i-will-file-a-bug-report.jpg)

It is a great help to the community if you report any bugs in the library
that you may find. We keep track of all open issues with the library
[here](https://github.com/dealii/dealii/issues). If you can, please try to
include a minimal failing example that can help us reproduce the problem.

## Making changes to the library
![fork on github](http://i.imgur.com/8jSBrnr.jpg)

To make a change to deal.II you should create a *fork* (through GitHub) of
the library and a separate *branch* (sometimes called a feature branch).

![merge into master](https://cdn.meme.am/instances/500x/37467752/does-not-simply-walk-into-mordor-boromir-one-does-not-simply-merge-into-master.jpg)

You can propose that your branch be combined with the rest of the library
by opening a *pull request*. This will give a chance for others to review
your code. 

![unit test](http://s2.quickmeme.com/img/03/0347c3efdc17cc1959d089f60b8b2fc267d9093caa8e8cb483bf476b58e63e45.jpg)

Each *pull request* should be accompanied by a small test that shows that 
the code works like intended. 

While this seems very formal, keeping all of the code review in one place
makes it easier to coordinate changes to the library (and there are usually
several people making changes to the library at once). This process is
described at length in
[lecture 32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html).
Please do not hesitate to ask questions about the workflow on the mailing
list if you are not sure what to do.

![coding style](https://cdn.meme.am/instances/500x/64488567/zoidberg-your-coding-style-is-bad-and-you-should-feel-bad.jpg)

Since deal.II is a fairly large project with lots of contributors we use a
set of
[coding conventions](https://www.dealii.org/developer/doxygen/deal.II/CodingConventions.html)
to keep the style of the source code in the library consistent. This
convention essentially consists of using `astyle` for indentation, camel
case for classes, and lower case letters with underscores for everything
else. If you are new to the project then we will work with you to ensure
your contributions are formatted with this style, so please do not think of
it as a road block if you would like to contribute some code.

The library is licensed under the GNU Lesser General Public License; while
you will retain copyright on your contributions, all changes to the core
library must be provided under this common license.
