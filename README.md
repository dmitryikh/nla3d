# nla3d

## Introduction

_nla3d_ (Non-Linear Analysis 3D) - is a program to solve mechanics non-linear finite element (FE)
problems. Or, better to say, it a set of procedures (C++ classes) that make it easy to implement new
FE, material model, constraint and other. It's written in C++ language massively using classes. Then
many low-level details are hidden from "end-user". This program is a good opportunity for students
and researchers to implement their own FE code by adding new element formulations, material models
to _nla3d_ framework.  The _nla3d_ code has been developing since 2011 year. First steps were made
by the author as a part of a course work. The name of the program was _nla_, because it was capable
to solve only 2D axisymmetric problems with Neo-Hookean hyperelastic material. Then the program was
improved for solving 3D problems in 2012 as a part of graduation work of the author. At those time
it was renamed to _nla3d_.

Currently, _nla3d_ has next features under the hood: 

  * It's written in c++, and therefore it has quite good performance. For example, on regular laptop
    an equilibrium iteration (single step of solution) for a FE model with ~ 100 000 degrees of
freedoms (DoFs) takes about 3-4 seconds.

  * For important parts of FE code _nla3d_ uses class hierarchy. For example, here is a base
    abstract class _Element_, and all particular FE formulations are derived from it. Therefore one
can easily integrate his own FE into _nla3d_.

  * Many commonly used procedures in FE methods are already implemented in _nla3d_ and they are
    hidden in classes realisations.

  * Here is two elements in _nla3d_:
    - *PLANE41* - 4-node element with u/p formulation (quasi linear
    interpolation of displacements and constant interpolation of hydrostatic pressure). This element
is suitable for model a behavior of model with nearly incompressible hyperelastic materials.
    - *SOLID81* the same as *PLANE41*, but for 3D models. All this elements support finite deformations with Total
  Lagrange approach.

  * It supports non-linear multi point constraints (MPC). _nla3d_ treats MPC as addition nonlinear
    equations to a global equation system. Here is rigid body MPC realization which can be useful in
a modeling of finite torsion.

  * Any hyperelastic materials described by invariants of a strain tensor can be used in the analysis.

  * Iterative approach to find an equilibrium state of a model for every load step. This is
    implemented by _Newton-Raphson method_ where non-linear global system of equations are linearized
and solved many times iteratively.

  * The mechanism of post-processors is implement. By this mechanism one can write his own code for
    results post-processing. For now, here are two already done post-processors: _VtkPostProcessor_ -
for writing results of solution in vtk file format. This format is widely used in scientific
visualization world.  One can analyze and plot results from vtk files in _Paraview_ open source
program.

What _nla3d_ can't do:

  * It's not all in all executable program. It more look like a set of classes or a library to help
    someone to implement his own FE code without needs to write every routine by himself.
Nevertheless, here is an executable _nla3d.exe_ that reads input data from Ansys's *.cdb file format
and solves the problem. _nla3d.exe_ also can be provided with many command line arguments to setup the
solver, boundary conditions, MPC, post-processors.

  * In current realization in _nla3d.exe_ it's impossible to use elements with different element types
    and with different material in a single model. But this limitation will be avoided in future
versions.

  * _nla3d_ isn't capable to solve transient problems.

## Requirements

_nla3d_ is written in C++ language. It uses C++11 features in the code. It is compiled ok with Visual
Studio 2013 compiler under Windows 7 and GNU g++ 4.7.2 under Debian. Actually, the Author didn't try
to compile it under different platforms and compilers. And here is non-zero probability that _nla3d_
will be compiled ok under another versions of compilers.

_nla3d_ uses CMake to deploy project files from the source code. CMake produces visual studio projects
under Windows and make files under linux-based OS. One can find more information about how to
compile _nla3d_ using CMake in the section *Compilation*.

_nla3d_ uses several external libraries:

  * Intel Math Kernel Library (MKL) to solve large scale systems of linear equations and to call BLAS
    matrix procedures (optional).
  * Eigen matrix library to provide rich and good optimized matrix operations.
  * Easylogging++ as powerfull and lightweight logging engine with the useful performance measurments
    feature.

## Compiling

To compile _nla3d_ first of all you need an environment and tools desribed in Requirements chapter.
Here is an example how to do it under linux-based OS. Firs of all we need to clone current project
and fetch all submodules:

```
git clone https://github.com/dmitryikh/nla3d
cd nla3d
git submodule update --init
```

Then you can deploy makefiles dy calling ```cmake``` program:

```
mkdir build
cd build
cmake ..
```

One will see the output like this:

```
...
-- Found MKL: /opt/intel  
-- Found EASYLOGGINGPP: /Users/foo/code/nla3d/site-src/easyloggingpp/src  
-- Found Eigen: /Users/foo/code/nla3d/site-src/eigen  
-- Configuring done
-- Generating done
...
```

That means that everything is ok with configuring.
If cmake wasn't successed to find were are some dependencies are, you need to fix manually. The
best way is to use cmake curses interface ```ccmake``` and provide missing pathes manually:

```
ccmake .
```

Then one can just launch the compilation:

```
cmake --build . --config Release 
```

And after launch the tests:

```
ctest . --C Release
```

If all test pass ok, in this case here is fully worked _nla3d_ binaries!

## Compiling options

Here are some _nla3d_ specific options that can affect on result binaries. This options are defined
in ```nla3d/CMakeLists.txt```:

```
set (nla3d_use_MKL OFF)
set (nla3d_multithreaded OFF)
set (nla3d_debug OFF)
set (nla3d_blas OFF)
```

_nla3d_ by default use ```math::GaussDenseEquationSolver``` to solve a system of linear equations. This is
very simple method named _Gaussian Elimination with partial pivoting_. It works fine for small
systems but when you need to solve some large FE system your should use other ```math::EquationSolver```.
The best candidate for now is ```math::PARDISO_equationSolver```. It uses _MKL_'s ```PARDISO(...)```
subroutines to solve large scale sparce linear equations systems. Hence, if you have _MKL_ dev
libraries installed on your side you can ```set (nla3d_use_MKL ON)``` to build with _MKL_ support.
Once you did this nla3d.exe will use ```math::PARDISO_equationSolver```.

Currently _nla3d_ internals doesn't benefit from multithread support. ```set (nla3d_multithreaded
ON)``` will only affect on _MKL_'s pardiso equation solver.

```set (nla3d_blas ON)``` will use blas subroutines to perform matrix calculations. Currently it's
not well tested, so it's recommended to keep this thing OFF.

```set (nla3d_debug ON)``` will only affect on compiler flags. See ```nla3d/CMakeLists.txt``` for
details.

 
## Writing new finite element

Here is just a try to help someone to get familiar with _nla3d_. In this chapter a pretty easy 3D
truss finite element is described and is incorporated into _nla3d_ library and then an excecutable
program is shown. 3D truss is a finite element with two spatial nodes connected with a line, this
line can be subjected only with compression or tension loads. The line has a few properties:
cross-section area, Young's module (stiffness of the material).

The implementation of a such FE one can find in ```src/lib/elements/TRUSS3.h``` and ```src/lib/elements/TRUSS3.cpp```.

Then ElementTRUSS3 is used to solve a very simple 2D problem described in ```src/main_truss.cpp```.

All this sources are massively commented to make it clear how to use it and how to implement new FE
into nla3d.

## Contacts

In case if you are interested in the project or if you have questions, please contact with by email:
khdmitryi ```at``` gmail.com
