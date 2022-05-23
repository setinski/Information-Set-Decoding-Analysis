# Information Set Decoding Analysis
In this project, we analyze the complexity of the classical and quantum information set decoding algorithms used to solve the generalized version of the syndrome decoding problem. For more details, refer to https://arxiv.org/abs/2104.12810.

## Project's dependencies

Before building and/or running the project, one needs to:  
1. [Download](https://www.boost.org/users/download/) **Boost** set of libraries and extract the archive: the path to the extracted folder will be your <BOOST_DIR> path.
2. [Download](https://www.mosek.com/downloads/) and [install](https://docs.mosek.com/9.2/cxxfusion/install-interface.html#testing-the-installation-and-compiling-examples) **MOSEK** software package: the path to the ../tools/.. subdirectory of **MOSEK** install path will be your <TOOLS_DIR> path.


## Build and Run the project

The project is built through the following steps:

1. Download the project. If downloaded as .zip archive, extract the directory from the downloaded archive: we denote this directory as <EXTRACTED_DIR>.
2. Navigate to the the subdirectory of downloaded (extracted) directory that corresponds to your operating system (namely, either Mac or Linux).
3. In the chosen subdirectory, find Makefile, open it in a text editor and edit the paths in Makefile so that they match the MOSEK and Boost paths on your machine (namely, <BOOST_DIR> and <TOOLS_DIR> from Makefile), as well as the platform (namely <PLATFORM> from Makefile).
4. Navigate to the downloaded (extracted) directory through the console and build the project using **make** directive.
  
For Linux:
<pre translate="no" dir="ltr" is-upgraded="">cd <EXTRACTED_DIR>/Linux
make
</pre>
For Mac:
<pre translate="no" dir="ltr" is-upgraded="">cd <EXTRACTED_DIR>/Mac
make
</pre>

If everything is done correctly, you will find the InformationSetDecodingAnalaysis executable in ./bin subdirectory of your downloaded (extracted) directory.

## Run the project

To run the project, navigate to the binary subdirectory and run the InformationSetDecoding executable:
<pre translate="no" dir="ltr" is-upgraded="">cd EXTRACTED_DIR/bin
./InformationSetDecoding
</pre>

If everything works correctly, through the console, you will be asked to choose the arguments of the syndrome decoding problem and the algorithm that solves it. In particular, you will choose between the Hamming and Lee weight and between Prange's, Dumer's and Wagner's algorithm, as well as if the selected algorithm runs in the classical or quantum regime. You will also choose the alphabet size of the observed problem.

In the resulting file (result.txt or result_quantum.txt), you will obtain the optimal values of the parameters of the selected algorithm when solving the hardest instance of the selected syndrome decoding problem, as well as running time of the chosen algorithm when the parameters are optimal.

### Authors
André Chailloux, INRIA Paris, COSMIQ project-team  
Thomas Debris-Alazard, INRIA Saclay, Grace project-team  
Simona Etinski, University of Paris, INRIA Paris, COSMIQ project-team  
