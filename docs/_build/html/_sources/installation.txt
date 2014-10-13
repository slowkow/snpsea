Installation
------------

**On Linux 64-bit, you may use the provided executable**

This runs on kernel 2.6.18 and newer:
https://github.com/slowkow/snpsea/releases

**Otherwise, you must build the executable from source**

The source code is available: https://github.com/slowkow/snpsea

**Mac:** To compile C++ code with the required dependencies, you need
XCode and MacPorts: http://guide.macports.org/#installing.xcode

Install the dependencies:

.. code:: bash

    # Ubuntu
    sudo apt-get install build-essential libopenmpi-dev libgsl0-dev

    # Mac
    #   First, install port (MacPorts): http://www.macports.org/
    #   Next, use it to install the dependencies:
    sudo port selfupdate && sudo port install gcc48 openmpi gsl

    # Broad Institute
    #   Add this line to ~/.my.bashrc or ~/.my.cshrc
    use .gcc-4.8.1 .openmpi-1.4 .gsl-1.14

Download and compile the code:

.. code:: bash

    #   Clone with git; easily get updates with 'git pull':
    git clone https://github.com/slowkow/snpsea.git
    cd snpsea

    #   If you don't have git:
    curl -LOk https://github.com/slowkow/snpsea/archive/master.zip
    unzip master.zip; cd snpsea-master

    cd src; make               #   Compile.
    cp ../bin/snpsea* ~/bin/   #   Copy the executables wherever you like.

C++ Libraries
~~~~~~~~~~~~~

To compile SNPsea, you will need a modern C++ compiler that supports
`c++0x <https://gcc.gnu.org/projects/cxx0x.html>`__ and the dependencies
listed below. I compiled successfully with gcc versions 4.6.3 (the default
version for Ubuntu 12.04) and 4.8.1.


`intervaltree <https://github.com/slowkow/intervaltree>`__

    a minimal C++ interval tree implementation

`Eigen <http://eigen.tuxfamily.org>`__

    Eigen is a C++ template library for linear algebra: matrices,
    vectors, numerical solvers, and related algorithms.

`OpenMPI <http://www.open-mpi.org>`__

    MPI is a standardized API typically used for parallel and/or
    distributed computing. Open MPI is an open source, freely available
    implementation.

`GSL - GNU Scientific Library <http://www.gnu.org/software/gsl>`__

    The GNU Scientific Library (GSL) is a numerical library for C and
    C++ programmers.

`GCC, the GNU Compiler <http://gcc.gnu.org>`__

    The GNU Compiler Collection is a compiler system produced by the GNU
    Project supporting various programming languages.

Python Packages
~~~~~~~~~~~~~~~

To plot visualizations of the results, you will need Python 2.7 and the
packages listed below.

**Instructions:** Install with `pip <http://www.pip-installer.org>`__:

::

    pip install docopt numpy pandas matplotlib

**Note:** The packages available on the Ubuntu repositories may be
outdated and might fail to work. So, avoid using ``apt-get`` for these
dependencies.

`docopt <http://docopt.org/>`__

    Command-line interface description language.

`numpy <http://www.numpy.org>`__

    NumPy is the fundamental package for scientific computing with
    Python.

`pandas <http://pandas.pydata.org>`__

    pandas is an open source, BSD-licensed library providing
    high-performance, easy-to-use data structures and data analysis
    tools for the Python programming language.

`matplotlib <http://matplotlib.org>`__

    matplotlib is a python 2D plotting library which produces
    publication quality figures in a variety of hardcopy formats and
    interactive environments across platforms.

**Note:** On a server with no display, please edit your
`matplotlibrc <http://matplotlib.org/users/customizing.html>`__ file to
use the ``Agg`` backend:

::

    perl -i -pe 's/^(\s*(backend).*)$/#$1\n$2:Agg/' ~/.matplotlib/matplotlibrc

Otherwise, you may see an error message like this:

::

    _tkinter.TclError: no display name and no $DISPLAY environment variable

R Packages
~~~~~~~~~~

Some visualizations use R and ggplot2 instead of Python and matplotlib.

**Instructions:** Start a session in R and run:

.. code:: r

    install.packages(c("data.table", "reshape2", "gap", "ggplot2"))


`data.table <http://cran.r-project.org/web/packages/data.table>`__

    Extension of data.frame for fast indexing, fast ordered joins, fast
    assignment, fast grouping and list columns.

`reshape2 <http://cran.r-project.org/web/packages/reshape2>`__

    Flexibly reshape data: a reboot of the reshape package.

`gap <http://cran.r-project.org/web/packages/gap>`__

    Genetic analysis package.

`ggplot2 <http://cran.r-project.org/web/packages/ggplot2>`__

    An implementation of the Grammar of Graphics.

