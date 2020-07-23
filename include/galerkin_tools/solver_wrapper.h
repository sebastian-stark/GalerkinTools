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


// TODO: internally UMFPACK interface of deal.II is used to invert the matrix related to the second block, replace this by something more efficient!
#ifdef DEAL_II_WITH_UMFPACK
#ifdef GALERKIN_TOOLS_WITH_UMFPACK

/**
 * A SolverWrapper for the direct solver from UMFPACK using directly the interface to UMFPACK (this allows for more flexibility of the solver)
 *
 * @note UMFPACK has its own license, independent of that of deal.II. If you
 * want to use the UMFPACK you have to accept that license.
 *
 * Code partially copied over from deal.II
 */
class BlockSolverWrapperUMFPACK2 : public SolverWrapper<dealii::Vector<double>, dealii::BlockVector<double>, TwoBlockMatrix<dealii::SparseMatrix<double>>, TwoBlockSparsityPattern>
{
private:

	/**
	 * Control array for the solver routines.
	 */
	mutable std::vector<double>
	control = std::vector<double>(UMFPACK_CONTROL);

	/**
	 * Info array for the solver routines.
	 */
	mutable std::vector<double>
	info = std::vector<double>(UMFPACK_INFO);

	/**
	* The size of the matrix A of the TwoBlockMatrix
	*/
	mutable unsigned int
	N;

	/**
	 * UMFPACK symbolic decomposition of A
	 */
	mutable void*
	symbolic_decomposition = nullptr;

	/**
	 * UMFPACK numeric decomposition of A
	 */
	mutable void*
	numeric_decomposition = nullptr;

	/**
	 * Array saying which row starts where in BlockSolverWrapperUMFPACK2::Ai
	 */
	mutable std::vector<SuiteSparse_long>
	Ap;

	/**
	 * Array storing the column indices of non-zero entries of A
	 */
	mutable std::vector<SuiteSparse_long>
	Ai;

	/**
	 * Array storing the numerical values of the non-zero entries of the matrix A
	 */
	mutable std::vector<double>
	Ax;

	/**
	 * Sets the data structures storing the matrix up (BlockSolverWrapperUMFPACK2::Ap, BlockSolverWrapperUMFPACK2::Ai, BlockSolverWrapperUMFPACK2::Ax).
	 *
	 * Also sets up the control parameters for the UMFPACK routines
	 *
	 * @param[in]	matrix		The matrix A
	 */
	void
	initialize_matrix(const SparseMatrix<double>& matrix)
	const;

	/**
	 * This function analysis the given matrix and sets up BlockSolverWrapperUMFPACK2::symbolic_decomposition.
	 * Whether this function does anything depends upon the value of BlockSolverWrapperUMFPACK2::analyze.
	 */
	void
	analyze_matrix()
	const;

	/**
	 * This function factorizes the matrix. The result is available from BlockSolverWrapperUMFPACK2::symbolic_decomposition.
	 */
	void
	factorize_matrix()
	const;

	/**
	 * Multiply the inverse of A by @p f and store the result into @p x. I.e., solve Ax=f for x.
	 * This uses the factors computed by BlockSolverWrapperUMFPACK2::factorize_matrix().
	 *
	 * @param[out]	x	solution
	 * @param[in]	f	rhs
	 */
	void
	vmult(	Vector<double>& 		x,
			const Vector<double>&	f)
	const;

public:

	/**
	 * Destructor.
	 */
	~BlockSolverWrapperUMFPACK2();

	/**
	 * Key indicating when to analyze the matrix structure:<br>
	 * 0 - before each factorization (default)<br>
	 * 1 - only before the next factorization (afterwards, BlockSolverWrapperUMFPACK2::analyze will be set to 2)<br>
	 * >=2 - do not recompute
	 *
	 * @warning		The user must ensure that BlockSolverWrapperUMFPACK2::analyze=0 or BlockSolverWrapperUMFPACK2::analyze=1
	 * 				whenever the sparsity pattern of the matrix has changed (or during the first call). Currently, no internal checking is performed.
	 */
	mutable unsigned int
	analyze = 0;

	/**
	 * level of diagnostic printing:<br>
	 * 0: no output<br>
	 * 1: error messages only<br>
	 * 2 or more: print status, whether or not an error occurred<br>
	 * 4 or more: also print the UMFPACK Copyright<br>
	 * 6 or more: also print the UMFPACK License
	 */
	unsigned int
	print_level = 0;

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

	DeclException2(	ExcUMFPACKError,
					std::string, int,
					<< "UMFPACK routine " << arg1 << " returned error status " << arg2 << ".");

};

#endif // GALERKIN_TOOLS_WITH_UMFPACK
#endif // DEAL_II_WITH_UMFPACK

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* INCREMENTALFE_SOLVERWRAPPER_H_ */
