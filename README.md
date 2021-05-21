# Information Set Decoding Analysis
In this project, we analyze the complexity of the classical and quantum information set decoding algorithms used to solve the generalized version of the syndrome decoding problem. For more details, refer to https://arxiv.org/abs/2104.12810.

## Requirements

Before building the project, one needs to have installed the following:  
1. [**Boost**](https://www.boost.org/users/download/) set of libraries;
2. [**MOSEK**](https://www.mosek.com/downloads/) software package;
3. [**CMake**](https://cmake.org/download/) software.

## Build the project

You can find a pre-built code of a project in the .\bin subdirectory of this project. However, it is built for Win64 in Release mode. For other platforms, it is necessary to build the project, and it can be done by performing the following steps:

1. download the project (if downloaded as .zip archive, extract the folder from the downloaded archive)
2. open the terminal (Linux)/command line(Windows)/(Mac) and navigate to the downloaded (extracted) folder
<pre translate="no" dir="ltr" is-upgraded="">cd PATH_TO_EXTRACTED_DIR
</pre>
3. build the project using CMake

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

If everything is done correctly, you will find ./bin subdirectory in your local directory.

## Run the project

To run the project, navigate to ./bin directory:
<pre translate="no" dir="ltr" is-upgraded="">cd PATH_TO_EXTRACTED_DIR/bin
</pre>

And run the InformationSetDecodingAnalysis.exe:
<pre translate="no" dir="ltr" is-upgraded="">InformationSetDecodingAnalysis.exe
</pre>

### Authors
André Chailloux, INRIA Paris, COSMIQ project-team
Thomas Debris-Alazard, INRRIA Saclay, Grace project-team
Simona Etinski, University of Paris, INRIA Paris, COSMIQ project-team
