The main purpose of the GalerkinTools library is to help with the
assembly of finite element systems for formulations obtained based on the Galerkin approach. Though the library must
of course make certain restrictions on the class of problems for which it can be applied, it is believed to
be useful for the solution of many practical problems especially in the field of solid mechanics. Special focus of
the library is on problems involving different domain portions (e.g. different materials) as well as on problems
involving unknown fields living on interfaces or boundaries.

The library currently requires deal.II 9.4, lapacke and numdiff to be installed on your system (and cmake must be able to find them).
To use the direct sequential solver interface of the library, deal.II must additionally be configured with UMFPACK.
To use the distributed parallel capabilities of the library, deal.II must be configured with MPI, PETSc (including MUMPS for the parallel solver) and p4est.

Installation of the library is through cmake:

1. place library source files into some folder /path/to/folder/GalerkinTools (you can use git clone https://github.com/sebastian-stark/GalerkinTools.git for this)
2. cd /path/to/folder/
3. mkdir build
4. cd build
5. cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ../GalerkinTools
6. make install
7. optionally set an environment variable to GALERKIN_TOOLS_DIR=/path/to/install/dir
8. optionally run the tests (first cd /path/to/folder/build, then ctest)

Acknowledgements:
The GalerkinTools library has been developed during a project supported by the Deutsche Forschungsgemeinschaft (DFG) under Grants STA 1593/1-1 and STA 1593/2-1.
