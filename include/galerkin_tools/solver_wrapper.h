// --------------------------------------------------------------------------
// Copyright (C) 2020 by Sebastian Stark
//
// This file is part of the GalerkinTools library
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef INCREMENTALFE_SOLVERWRAPPER_H_
#define INCREMENTALFE_SOLVERWRAPPER_H_

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <galerkin_tools/config.h>
#include <galerkin_tools/two_block_sparsity_pattern.h>
#include <galerkin_tools/two_block_matrix.h>

#ifdef GALERKIN_TOOLS_WITH_UMFPACK
#include <umfpack.h>
#endif // GALERKIN_TOOLS_WITH_UMFPACK

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * This class essentially wraps some functions related to the solver used for the problem
 * and the corresponding data types into a common interface.
 *
 * The SolverWrapper class inherits from Subscriptor in order to be
 * able to check that SolverWrapper objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam		SolutionVectorType		The type used for the solution vector
 *
 * @tparam		RHSVectorType			The type used for the right hand side
 *
 * @tparam		MatrixType				The type used for the system matrix
 *
 * @tparam		SparsityPatternType		The sparsity pattern used for initialization of the matrix
 */

template<class SolutionVectorType, class RHSVectorType, class MatrixType, class SparsityPatternType>
class SolverWrapper : public dealii::Subscriptor
{
public:

	/**
	 * The solve function for the system of linear equations in the stretched form
	 *
	 * @param[in]	K_stretched		The stretched system_matrix
	 *
	 * @param[out]	solution		The stretched solution vector
	 *
	 * @param[in]	f_stretched		The stretched right hand side vector
	 *
	 * @param[in]	symmetric		A bool indicating whether @p system_matrix is symmetric
	 */
	virtual
	void
	solve(	const MatrixType&		K_stretched,
			SolutionVectorType&		solution,
			const RHSVectorType&	f_stretched,
			const bool				symmetric)
	const = 0;

	/**
	 * The destructor of SolverWrapper essentially checks before destruction that the
	 * SolverWrapper object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual
	~SolverWrapper();
};

/**
 * A SolverWrapper for the direct solver from UMFPACK using the UMFPACK interface of deal.II
 */
class SolverWrapperUMFPACK : public SolverWrapper<dealii::Vector<double>, dealii::Vector<double>, dealii::SparseMatrix<double>, SparsityPattern>
{
public:

	/**
	 * @copydoc SolverWrapper::solve
	 */
	virtual
	void
	solve(	const dealii::SparseMatrix<double>&	K_stretched,
			dealii::Vector<double>&				solution,
			const  dealii::Vector<double>&		f_stretched,
			const bool							symmetric = false)
	const;

};

/**
 * A SolverWrapper for the direct solver from UMFPACK for the case of a TwoBlockMatrix using the UMFPACK interface of deal.II
 */
class BlockSolverWrapperUMFPACK : public SolverWrapper<dealii::Vector<double>, dealii::BlockVector<double>, TwoBlockMatrix<dealii::SparseMatrix<double>>, TwoBlockSparsityPattern>
{
public:

	/**
	 * @copydoc SolverWrapper::solve
	 */
	virtual
	void
	solve(	const TwoBlockMatrix<dealii::SparseMatrix<double>>&	K_stretched,
			dealii::Vector<double>&								solution,
			const  dealii::BlockVector<double>&					f_stretched,
			const bool											symmetric = false)
	const;

};

#ifdef DEAL_II_WITH_PETSC
#ifdef DEAL_II_PETSC_WITH_MUMPS
#ifdef DEAL_II_WITH_MPI

/**
 * A SolverWrapper for the PETSc MUMPS solver, which works for parallel computations.
 */
class SolverWrapperPETSc : public SolverWrapper<dealii::LinearAlgebra::distributed::Vector<double>, dealii::PETScWrappers::MPI::BlockVector, parallel::TwoBlockMatrix<dealii::PETScWrappers::MPI::SparseMatrix>, TwoBlockSparsityPattern>
{
public:

	/**
	 * @copydoc SolverWrapper::solve
	 */
	virtual
	void
	solve(	const parallel::TwoBlockMatrix<dealii::PETScWrappers::MPI::SparseMatrix>&	K_stretched,
			dealii::LinearAlgebra::distributed::Vector<double>&							solution,
			const dealii::PETScWrappers::MPI::BlockVector&								f_stretched,
			const bool																	symmetric = false)
	const;
};

#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_PETSC_WITH_MUMPS
#endif // DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_PETSC
#ifdef DEAL_II_WITH_MPI

/**
 * A SolverWrapper for the PETSc MUMPS solver, which works for parallel computations.
 */
class SolverWrapperPETScIterative : public SolverWrapper<dealii::LinearAlgebra::distributed::Vector<double>, dealii::PETScWrappers::MPI::BlockVector, parallel::TwoBlockMatrix<dealii::PETScWrappers::MPI::SparseMatrix>, TwoBlockSparsityPattern>
{
public:

	/**
	 * @copydoc SolverWrapper::solve
	 */
	virtual
	void
	solve(	const parallel::TwoBlockMatrix<dealii::PETScWrappers::MPI::SparseMatrix>&	K_stretched,
			dealii::LinearAlgebra::distributed::Vector<double>&							solution,
			const dealii::PETScWrappers::MPI::BlockVector&								f_stretched,
			const bool																	symmetric = false)
	const;
};

#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_WITH_MPI

#ifdef GALERKIN_TOOLS_WITH_PARDISO

/** PARDISO prototypes */
extern "C" void pardisoinit(void*, int*, int*, int*, double*, int *);
extern "C" void pardiso(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, double*, double*, int*, double*);
extern "C" void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
extern "C" void pardiso_chkvec(int*, int*, double*, int*);
extern "C" void pardiso_printstats(int*, int*, double*, int*, int*, int*, double*, int*);

/**
 * A SolverWrapper for the PARDISO solver included in the Intel MKL library.
 */
class BlockSolverWrapperPARDISO : public SolverWrapper<dealii::Vector<double>, dealii::BlockVector<double>, TwoBlockMatrix<dealii::SparseMatrix<double>>, TwoBlockSparsityPattern>
{
public:

	/**
	 * @copydoc SolverWrapper::solve
	 */
	virtual
	void
	solve(	const TwoBlockMatrix<dealii::SparseMatrix<double>>&	K_stretched,
			dealii::Vector<double>&								solution,
			const  dealii::BlockVector<double>&					f_stretched,
			const bool											symmetric = false)
	const;
};

#endif // GALERKIN_TOOLS_WITH_PARDISO


GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* INCREMENTALFE_SOLVERWRAPPER_H_ */
