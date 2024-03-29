# GiNaCDE- an NLPDE solver
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03885/status.svg)](https://doi.org/10.21105/joss.03885)


GiNaCDE is an NLPDE (Nonlinear Partial Differential Equation) and NLODE (Nonlinear Ordinary Differential Equation) solver written in C++. It has three different methods: F-expansion, modified F-expansion (mF-expansion), and first integral method (FIM) for solving NLPDEs and NLODEs. It can be used to get exact analytical traveling-wave solutions of a wide variety of NLPDEs arising in different scientific community fields.
The library has been designed to solve the NLPDEs, which have the following general form

![](img/Eqn.png)

where ![](img/Eqn1.png) are independent variables, ![](img/Eqn2.png) is dependent variable, ![](img/Eqn3.png) are the parameters. Here *F* must be a polynomial about *u* and its derivatives.
GiNaCDE always transforms the NLPDE into an NLODE with respect to the traveling-wave coordinate ![](img/Eqn4.png) using a traveling-wave transformation.
F-expansion and modified F-expansion methods can be applied to higher-order NLPDEs. But, FIM is applicable to an NLPDE when its transformed NLODE with respect to the traveling-wave coordinate ![](img/Eqn5.png) is second-order only.
However, there is no guarantee that the library always gives the complete solutions of all NLPDEs of the above form. Sometimes, the library may fail to give solutions due to the complexity of the problems.

**N.B.: Now, GiNaCDE (>=v1.6.0) determines the solutions of differential equations without assuming all the constant parameters are strictly real and positive, i.e., all the constant parameters can be real and positive or real and negative numbers.**

## Features
Some interesting features of GiNaCDE are

  * Available solution methods: F-expansion method, modified F-expansion method, and first integral method.
  * It can solve  NLPDEs or NLODEs that contain complex functions.
  * It can tackle the non-polynomial form of *h (X)* in the case of FIM.
  * It can integrate an integrable NLPDE or NLODE when possible, and the generated integration constants can be assigned with the values in the user choice.
  * For differential equations with parameters (![](img/Eqn3.png)), it can determine the conditions on the parameters to obtain exact solutions.
  * The exact analytical solutions of NLPDEs or NLODES with calculating steps are saved in a text file written in `MAPLE`, `MATHEMATICA`, or `GiNaC` language.
  * It has a friendly Graphical User Interface (GUI).

## Example usage  
This is a brief example that computes exact solutions of the following one-dimensional cubic nonlinear Schrödinger (NLS) equation:

![](img/Eqn6.png)

where *p* and *q* are are non-zero real constants and *u(x,t)* is a complex-valued function depends on the variables *t,x*.
For a more detailed introduction, please refer to the [documentation](doc/documentation.pdf) file.  
```c++
/** @file NLS.cpp
 *
 *   This program solves one dimensional cubic nonlinear Schr\"odinger (NLS) equation:
        Iu_t-pu_{xx}+q{|u|}^2u=0, 
 **/

#include <GiNaCDE/GiNaCDE.h>

int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), k_0=reader("k_0"), k_1=reader("k_1"),
             p_0=reader("p_0"), p_1=reader("p_1"), A_0=reader("A_0"),A_2=reader("A_2"),
             p=reader("p"),q=reader("q");

    const ex pde = I*Diff(u,t,1) - p*Diff(u,x,2) + q*u*u*conjugate(u);

    depend(u, {t, x});

    output=maple;
    twcPhase=lst{lst{k_0,k_1},lst{p_0,p_1}};
    degAcoeff=lst{2,A_0,0,A_2};
    ASolve=false;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="NLS_Fexp(maple).txt";
    desolve(pde,{u},F_expansion);
    output=ginac;
    filename="NLS_Fexp_ginac.txt";
    desolve(pde,{u},F_expansion);

    output=mathematica;
    filename="NLS_FIM.txt";
    desolve(pde, {u}, FIM);

    return 0;

}
```
After compiling and running the above program, exact solutions with calculating steps are saved in the text files [NLS_Fexp(maple).txt](examples/NLS_Fexp(maple).txt), [NLS_Fexp_ginac.txt](examples/NLS_Fexp_ginac.txt) and [NLS_FIM.txt](examples/NLS_FIM.txt). 

## Installation

 
### External dependencies
GiNaCDE V1.6.0 requires the packages [CLN >= 1.3.4](http://www.ginac.de/CLN/), [GiNaC >= 1.8.1](https://www.ginac.de/archives/) and [GTK+ 3.xx](https://download-fallback.gnome.org/sources/gtk+/3.24/) (this library is optional and is used to build the GUI version of the GiNaCDE library). 

##### For Linux/MacOS machines:
All the dependencies are available via the most common package managers `APT` on Ubuntu or Debian. Additionally, all dependencies can also be retrieved on `macOS` via the most common package managers [Homebrew](https://brew.sh/) and [MacPorts](https://ports.macports.org/). For example, CLN, GiNaC and GTK3 are installed via APT through

```apt install libcln-dev libginac-dev libgtk-3-dev```

and on Homebrew via

```brew install cln ginac gtk+3```

##### For Windows machines:

 We have to install the libraries  [CLN >= 1.3.4](http://www.ginac.de/CLN/), [GiNaC >= 1.8.1](https://www.ginac.de/archives/) using the following commands:  
```        
        $ ./configure
        $ make
        $ make install
```
 The library  [GTK+ 3.xx](https://download-fallback.gnome.org/sources/gtk+/3.24/) can be easily installed using [MSYS2](https://msys2.github.io/), which provides a UNIX-like environment for Windows. It provides packages for many software applications and libraries, including the GTK stack. To install GTK3 and its dependencies, open a MSYS2 shell, and run:
```
        $ pacman -S mingw-w64-x86_64-gtk3
```   

### Compiling and installing GiNaCDE 
Compilations are done using the tools [CMake](https://cmake.org/download/) `>= 3.1`, pkg-config `>=0.29.2`, and of course, a compiler that supports `>= C++11`.
We suggest to use the C++ compiler from the GNU compiler collection, `GCC >= 4.9`. GiNaCDE can then be compiled using the commands: 
```
     $ mkdir build-dir # generate a separate directory
     $ cd build-dir
     $ cmake -DGINACDE_GUI_BUILD=on <path-to-source> # generate Makefiles
     $ make
     $ make install
```
To install the software locally (e.g. into `~/.local` or similar), one needs to add `-DCMAKE_INSTALL_PREFIX:PATH=~/.local` to the `cmake` call, i.e.
```
cmake -DGINACDE_GUI_BUILD=on -DCMAKE_INSTALL_PREFIX:PATH=~/.local <path-to-source> # generate Makefiles
```

A successful compilation will lead to the creation of libraries, executables of gtools, and GiNaCDE-GUI.`gtools` is a console application. GiNaCDE-GUI is the graphical user interface (GUI) of GiNaCDE library. If you do not want to build GiNaCDE-GUI, use the following option:
```
    -DGINACDE_GUI_BUILD=off
```
We have checked that the source files are successfully compiled on Windows
platform using MSYS2 (https://www.msys2.org) when we use the command

	cmake -G "MSYS Makefiles" -DGINACDE_GUI_BUILD=on <path-to-source>
#### Automated Tests
The [`test`](test/) folder contains tests. To run the tests, execute one of the following commands:
```
    $ make test
```
or 
```
    $ ctest
```
under `build-dir` created earlier for building GiNaCDE.
These automated tests verify the functionalities of the software. 

#### Precompiled GiNaCDE-GUI for Windows
We have provided a pe-compiled GiNaCDE-GUI, which can be downloaded from [here](https://sourceforge.net/projects/ginacde). The GiNaCDE-GUI has been compiled on Windows 10 OS using [`MSYS2`](https://www.msys2.org), GCC 10.3.0, GTK+ 3.24.30, CLN 1.3.6 and GiNaC 1.8.1. The precompiled software is compatible with 32-bit and 64-bit Windows 10 OS.



## Execution
GiNaCDE library can be executed in C++ code with GNU compiler collection, `GCC >= 4.9`. 
To run the GiNaCDE library from the GCC compiler, use the following command:
```  
    $ g++ -std=c++11 -Wall -g example.cpp -o example -lcln -lginac -lGiNaCDE
```   
To run `GiNaCDE-GUI`, `gtools` just click on `GiNaCDE_gui.exe`, `gtools.exe` files respectively. Then GiNaCDE-GUI is executed in a GUI framework, 
but `gtools` is executed in a console.  

If we want to create a new CMake project that uses GiNaCDE, we need to link the GiNaCDE library against the project executables. GiNaCDE provides a `pkg-config` configuration. So we currently need to do the following:

```
find_package(PkgConfig)
pkg_search_module(GiNaCDE REQUIRED IMPORTED_TARGET GiNaCDE>=1.0)
add_executable(my_example my_example.cpp)
target_link_libraries(my_example PRIVATE GiNaCDE)
```

### Output
GiNaCDE prints all the output results in a separate text ('.txt') file.
Besides this, the solutions of the NLPDE are collected by the programming variables *solutionClt* and *constraints*.

### Checking the solutions of Diff. Equ.
We can easily verify the solutions returned by GiNaCDE by substituting the solutions back into the differential equation. Following this substitution method, the solutions of Differential Equations given in the text file derived by GiNaCDE, can be easily checked by the software: Maple, Mathematica, and GiNaCDE. To illustrate the procedures for checking the solutions, we have provided some output text files [`NLS_Fexp(maple).txt`](examples/NLS_Fexp(maple).txt), [`NLS_Fexp(ginac).txt`](examples/NLS_Fexp(ginac).txt), [`KDV_FIM2.txt`](examples/KDV_FIM2.txt), [`gardner_Fex.txt`](examples/gardner_Fex.txt), [`cahnAllen_mF.txt`](examples/cahnAllen_mF.txt), [`Painlev_FIMextravar.txt`](examples/Painlev_FIMextravar.txt), and the corresponding checking files [`checkSolu_NLS_Fexp(maple).mw`](examples/checkSolu_NLS_Fexp(maple).mw), [`checkSolu_NLS_Fexp(ginac).cpp`](examples/checkSolu_NLS_Fexp(ginac).cpp), [`checkSolu_KDV_FIM2.mw`](examples/checkSolu_KDV_FIM2.mw), [`checkSolu_gardner_Fex.nb`](examples/checkSolu_gardner_Fex.nb), [`checkSolu_cahnAllen_mF.cpp`](examples/checkSolu_cahnAllen_mF.cpp), [`checkSolu_Painlev_FIMextravar.cpp`](examples/checkSolu_Painlev_FIMextravar.cpp) which explain how to test the solutions using Maple (Maple 2019), Mathematica (Mathematica 9), and GiNaCDE software.

**Caution:** Currently, GiNacDE is unable to check all the solutions reported by GiNaCDE due to some simplification problems. I hope this problem can be fixed in the future release of GiNaCDE.
Now to verify the solutions, I recommend to use Maple or Mathematica software. 


### Additional notes
We should note that we can obtain different results in output files after each running session of a GiNaCDE program. 
This happens because of the GiNaC library.
Because GiNaC assigns a unique (hidden) serial number for each newly created symbol object and GiNaC uses this unique serial number instead of its name for algebraic manipulations. The serial number for the same name of the symbol may be changed in each running session of the GiNaC program. As a result, the symbols in the same algebraic expressions may be ordered differently during each running session of the GiNaC program. This happens because to order the symbols of an algebraic expression GiNaC internally uses a comparison predicate, called *ex_is_less*, which uses an internal symbol id counter. 


## Examples
The [`examples`](examples/) folder contains all the examples which solve some NLPDEs, such as, Eckhaus equation, Seventh-order Sawada-Kotara equations, Fifth-order Generalized Korteweg–De Vries (KdV) equation, Perturbed nonlinear Schrödinger (NLS) Equation with Kerr Law Nonlinearity, Kudryashov-Sinelshchikov Equation, etc.
 
To compile the examples, move to the `build-dir` created earlier for building GiNaCDE, and execute
```
$ make examples
```
The executables will be placed into the `build-dir/bin` directory.


## Documentation: 
The documentation for GiNaCDE is available [`here`](doc/documentation.pdf).
The short tutorials on `GiNaCDE-GUI` and `gtools` are also available [`here`](doc/GiNaCDE_guiTutorial.pdf) and [`here`](doc/gtoolsTutorial.pdf), respectively.
    


    
## Contributions and bug reports
Contributions to this project are very welcome.
If you wish to contribute a new feature, you can do this by forking the GiNaCDE repo and creating a branch. Apply your code changes to the branch on your fork. When you're done, submit a [pull request](https://github.com/mithun218/GiNaCDE/pulls) to merge your fork into master branch with a tag "enhancement", and the proposed changes can be discussed there. 

If you encounter a bug, please open a new [issue](https://github.com/mithun218/GiNaCDE/issues/new) on the GitHub repository to report the bug, and tag it "bug".
Please provide sufficient information to reproduce the bug and include as much information as possible that can be helpful for fixing it.

## Citation
If you use this software in your research works, please cite the following article:

Bairagi, M., (2022). GiNaCDE: the high-performance F-expansion and First Integral Methods with C++ library for solving Nonlinear Differential Equations. Journal of Open Source Software, 7(72), 3885, https://doi.org/10.21105/joss.03885