# nla3d ## Introduction

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

  * Here is two elements in _nla3d_: *PLANE41* - 4-node element with u/p formulation (quasi linear
    interpolation of displacements and constant interpolation of hydrostatic pressure). This element
is suitable for model a behavior of model with nearly incompressible hyperelastic materials. *SOLID81*
- the same as *PLANE41*, but for 3D models. All this elements support finite deformations with Total
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

