# Information Set Decoding Analysis
In this project, we analyze the complexity of the classical and quantum information set decoding algorithms used to solve the generalized version of the syndrome decoding problem. For more details, refer to https://arxiv.org/abs/2104.12810.

## Project's dependencies

Before building and/or running the project, one needs to:  
1. [Download](https://www.boost.org/users/download/) **Boost** set of libraries and extract the archive:  
The path to your extracted folder will be your BOOST_DIR path.
2. [Download](https://www.mosek.com/downloads/) and [install](https://docs.mosek.com/9.2/cxxfusion/install-interface.html#testing-the-installation-and-compiling-examples) **MOSEK** software package:  
During installation, you will be only concerned with the section entitled "Manual compilation (all platforms)". Other sections are covered by the CMake.txt available in ./src directory of this repository.    
3. [Download](https://cmake.org/download/) and install **CMake** software.


## Build the project

You can find a pre-built code of a project in the ./bin/Windows/Release subdirectory of this project. However, it is built for Win64 in Release mode. For other platforms, it is necessary to build the project, and it can be done by performing the following steps:

1. download the project (if downloaded as .zip archive, extract the folder from the downloaded archive),
2. navigate to the downloaded (extracted folder), open CMake.txt in a text editor and edit the paths in CMake.txt so that they match the MOSEK and Boost paths in your system (namely, BOOST_DIR, MOSEK_HEADER_DIR, MOSEK_LIB, and FUSION_LIB),
3. open Command Promt (Windows)/Terminal (Linux and Mac) and navigate to the downloaded (extracted) folder through the console:
<pre translate="no" dir="ltr" is-upgraded="">cd PATH_TO_EXTRACTED_DIR
</pre>
4. build the project using CMake:

 - create **build** directory and navigate to it:
<pre translate="no" dir="ltr" is-upgraded="">mkdir build
cd build
</pre>

- run CMake to configure the project and generate a native build system:
<pre translate="no" dir="ltr" is-upgraded="">cmake ../src  
</pre>

- call the build system to build the project:
<pre translate="no" dir="ltr" is-upgraded="">cmake --build . --config Release
</pre>

If everything is done correctly, you will find the InformationSetDecodingAnalaysis.exe in ./bin/Release subdirectory of your local directory.

## Run the project

To run the project:  
1. open the consol application and navigate to the directory containing executable (the executable should be stored in PATH_TO_EXTRACTED_DIR/bin/Release):
<pre translate="no" dir="ltr" is-upgraded="">cd PATH_TO_EXE_DIR
</pre>

2. run the InformationSetDecodingAnalysis.exe:
<pre translate="no" dir="ltr" is-upgraded="">InformationSetDecodingAnalysis.exe
</pre>

If everything works correctly, through the console, you will be asked to choose the arguments of the syndrome decoding problem and the algorithm that solves it. In particular, you will choose between the Hamming and Lee weight and between Prange's, Dumer's and Wagner's algorithm, as well as if the selected algorithm runs in the classical or quantum regime. You will also choose the alphabet size of the observed problem. In the resulting file (result.txt or result_quantum.txt), you will obtain the optimal values of the parameters of the selected algorithm when solving the hardest instance of the selected syndrome decoding problem. 

### Authors
André Chailloux, INRIA Paris, COSMIQ project-team  
Thomas Debris-Alazard, INRIA Saclay, Grace project-team  
Simona Etinski, University of Paris, INRIA Paris, COSMIQ project-team  
