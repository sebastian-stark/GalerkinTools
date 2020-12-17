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

#include <umfpack.h>

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
	= 0;

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
			const bool							symmetric = false);

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
			const bool											symmetric = false);

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
			const bool																	symmetric = false);
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
			const bool																	symmetric = false);
};

#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_WITH_MPI

// TODO: internally UMFPACK interface of deal.II is used to invert the matrix related to the second block, replace this by something more efficient!
#ifdef GALERKIN_TOOLS_WITH_PARDISO
#ifdef DEAL_II_WITH_UMFPACK

/** PARDISO prototypes */
extern "C" void pardisoinit(void*, int*, int*, int*, double*, int *);
extern "C" void pardiso(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, double*, double*, int*, double*);
extern "C" void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
extern "C" void pardiso_chkvec(int*, int*, double*, int*);
extern "C" void pardiso_printstats(int*, int*, double*, int*, int*, int*, double*, int*);
extern "C" void pardiso_residual(int*, int*, double*, int*, int*, double*, double*, double*, double*, double*);

/**
 * A SolverWrapper for the PARDISO solver.
 *
 * @note PARDISO has its own license, independent of that of GalerkinTools. If you
 * want to use PARDISO you have to accept that license and obtain your own copy of PARDISO.
 *
 * Code partially copied over from deal.II
 */
class BlockSolverWrapperPARDISO : public SolverWrapper<dealii::Vector<double>, dealii::BlockVector<double>, TwoBlockMatrix<dealii::SparseMatrix<double>>, TwoBlockSparsityPattern>
{
private:

	/**
	 * stores whether PARDISO has been initialized.
	 */
	bool
	initialized = false;

	/**
	 * PARDISO control parameters -> refer to PARDISO user manual for details
	 */
	int
	iparm[64];

	/**
	 * PARDISO control parameters for iterative solvers -> here only a dummy
	 */
	double
	dparm[64];

	/**
	 * Pointers to internal solver memory
	 */
	void*
	pt[64];

	/**
	* The size of the matrix A of the TwoBlockMatrix
	*/
	int
	N;

	/**
	 * Array saying which row starts where in BlockSolverWrapperPARDISO::Ai
	 */
	std::vector<int>
	Ap;

	/**
	 * Array storing the column indices of non-zero entries of A
	 */
	std::vector<int>
	Ai;

	/**
	 * Array storing the numerical values of the non-zero entries of the matrix A
	 */
	std::vector<double>
	Ax;

	/**
	 *  Array storing the column scaling matrix if apply_scaling==2
	 */
	dealii::Vector<double>
	C_inv;

	/**
	 *  Array storing the row scaling matrix if apply_scaling==2
	 */
	dealii::Vector<double>
	R_inv;

	/**
	 *  Maximum number of numerical factorizations to keep in memory (PARDISO internal).
	 */
	int
    maxfct = 1;

	/**
	 * Number of right hand sides to solve for in one step
	 */
	int
	nrhs = 1;

	/**
	 *  Matrix number (PARDISO internal).
	 */
	int
	mnum = 1;

	/**
	 * diagnostic printing (PARDISO internal):<br>
	 * 0: no output<br>
	 * 1: print output
	 */
    int
	msglvl;

    /**
     * Vector for residual calculation
     */
    Vector<double>
    res;

	/**
	 * Sets the data structures storing the matrix up (BlockSolverWrapperPARDISO::Ap, BlockSolverWrapperPARDISO::Ai, BlockSolverWrapperPARDISO::Ax).
	 *
	 * Also sets up the control parameters for the PARDISO routines
	 *
	 * @param[in]	matrix		The matrix A
	 */
	void
	initialize_matrix(	const SparseMatrix<double>& matrix);

	/**
	 * This function analyses the given matrix.
	 * Whether this function does anything depends upon the value of BlockSolverWrapperPARDISO::analyze.
	 */
	void
	analyze_matrix();

	/**
	 * This function factorizes the matrix.
	 */
	void
	factorize_matrix();

	/**
	 * Multiply the inverse of A by @p f and store the result into @p x. I.e., solve Ax=f for x.
	 * This uses the factors computed by BlockSolverWrapperPARDISO::factorize_matrix().
	 *
	 * @param[out]	x	solution
	 * @param[in]	f	rhs
	 */
	void
	vmult(	Vector<double>& 		x,
			const Vector<double>&	f);

	/**
	 * @return	Return the matrix type (in PARDISO convention - the return is not equal to BlockSolverWrapperPARDISO::matrix_type)
	 */
	int get_matrix_type();

	/**
	 * Scales the system and computes BlockSolverWrapperPARDISO::C_inv and BlockSolverWrapperPARDISO::f_scaled
	 *
	 * @param[in]	K_stretched		the matrix
	 *
	 * @param[in]	f_stretched		the r.h.s.
	 */
	void
	scale_system(	TwoBlockMatrix<dealii::SparseMatrix<double>>*	K_stretched,
					dealii::BlockVector<double>*					f_stretched);

	/**
	 * Scales the solution by multiplying it by BlockSolverWrapperPARDISO::C_inv
	 *
	 * and undoes the scaling of K and f
	 *
	 * @param[inout]	solution	the solution vector
	 *
	 * @param[inout]	K_stretched	the original matrix
	 *
	 * @param[inout]	f_stretched	the original r.h.s.
	 */
	void
	scale_solution(	dealii::Vector<double>&							solution,
					TwoBlockMatrix<dealii::SparseMatrix<double>>*	K_stretched,
					dealii::BlockVector<double>*					f_stretched)
	const;

public:

	/**
	 * Destructor.
	 */
	~BlockSolverWrapperPARDISO();

	/**
	 * Key indicating when to analyze the matrix structure:<br>
	 * 0 - before each factorization (default)<br>
	 * 1 - only before the next factorization (afterwards, BlockSolverWrapperPARDISO::analyze will be set to 2)<br>
	 * >=2 - do not recompute
	 *
	 * @warning		The user must ensure that BlockSolverWrapperPARDISO::analyze=0 or BlockSolverWrapperPARDISO::analyze=1
	 * 				whenever the sparsity pattern of the matrix has changed (or during the first call). Currently, no internal checking is performed.
	 *
	 * @warning		In case that scaling==1, the scaling of the matrix is only recomputed if BlockSolverWrapperPARDISO::analyze <2. Otherwise, the most recently computed scaling will be used.
	 */
	unsigned int
	analyze = 0;

	/**
	 * diagnostic printing:<br>
	 * 0: no output<br>
	 * 1: print output
	 */
	unsigned int
	print_level = 0;

	/**
	 * type of the matrix:<br>
	 * 0: unsymmetric<br>
	 * 1: structurally symmetric<br>
	 * 2: symmetric indefinite<br>
	 * >=3: symmetric positive definite
	 */
	unsigned int
	matrix_type = 0;

	/**
	 * ordering method:<br>
	 * 0: AMD<br>
	 * 2: nested dissection (METIS)
	 */
	unsigned int
	ordering_method = 2;

	/**
	 * whether to compute a scaling of the matrix during analysis:<br>
	 * 0: no scaling<br>
	 * 1: matched scaling of PARDISO package<br>
	 * 2: iterative scaling method of D. Ruiz (D.Ruiz: A Scaling Algorithm to Equilibrate Both Rows and Columns Norms in Matrices, 2001)
	 */
	unsigned int
	apply_scaling = 0;

	/**
	 * relative tolerance for the iterative procedure if apply_scaling == 2 (1 is the highest tolerance and 0 the lowest tolerance)
	 *
	 * For 0.0 only a single iteration is performed
	 */
	double
	iterative_scaling_tolerance = 0.9;

	/**
	 * pivoting type for symmetric indefinite matrices<br>
	 * 0: 1x1<br>
	 * 1: 2x2
	 */
	unsigned int
	pivoting_method = 1;

	/**
	 * maximum number of iterative refinements
	 */
	unsigned int
	n_iterative_refinements = std::numeric_limits<int>::max();

	/**
	 * value of pivot to be added to tiny pivots (if e.g. 10, than 1e-10 is added to tiny pivots)
	 */
	unsigned int
	pivot_perturbation = 8;

	/**
	 * If @p true: use user defined permutation
	 */
	bool
	user_perm = false;

	/**
	 * Maximum allowable residual
	 */
	double
	res_max = 1e-8;

	/**
	 * indicates whether to use default parameters of PARDISO. If @p true, BlockSolverWrapperPARDISO::ordering_method, BlockSolverWrapperPARDISO::apply_scaling, BlockSolverWrapperPARDISO::pivoting_method,
	 * BlockSolverWrapperPARDISO::n_iterative_refinements, BlockSolverWrapperPARDISO::pivot_perturbation, BlockSolverWrapperPARDISO::user_perm are ignored.
	 */
	bool
	use_defaults = true;

	/**
	 * user defined permutation vector
	 */
	std::vector<int>*
	perm;

	/**
	 * @copydoc SolverWrapper::solve
	 */
	virtual
	void
	solve(	const TwoBlockMatrix<dealii::SparseMatrix<double>>&	K_stretched,
			dealii::Vector<double>&								solution,
			const  dealii::BlockVector<double>&					f_stretched,
			const bool											symmetric = false);

	DeclException2(	ExcPARDISOError,
					std::string, int,
					<< "PARDISO routine " << arg1 << " returned error status " << arg2 << ".");

	DeclException1(	ExcPARDISORes,
					double,
					<< "PARDISO residual too large: res = " << arg1 << ".");

};

#endif // DEAL_II_WITH_UMFPACK
#endif // GALERKIN_TOOLS_WITH_PARDISO

// TODO: internally UMFPACK interface of deal.II is used to invert the matrix related to the second block, replace this by something more efficient!
#ifdef DEAL_II_WITH_UMFPACK

/**
 * A SolverWrapper for the direct solver from UMFPACK using directly the interface to UMFPACK (this allows for more flexibility of the solver compared to BlockSolverWrapperUMFPACK, which uses the deal.II
 * interface to UMFPACK)
 *
 * @note UMFPACK has its own license, independent of that of GalerkinTools. If you
 * want to use UMFPACK you have to accept that license.
 *
 * Code partially copied over from deal.II
 */
class BlockSolverWrapperUMFPACK2 : public SolverWrapper<dealii::Vector<double>, dealii::BlockVector<double>, TwoBlockMatrix<dealii::SparseMatrix<double>>, TwoBlockSparsityPattern>
{
private:

	/**
	 * Control array for the solver routines.
	 */
	std::vector<double>
	control = std::vector<double>(UMFPACK_CONTROL);

	/**
	 * Info array for the solver routines.
	 */
	std::vector<double>
	info = std::vector<double>(UMFPACK_INFO);

	/**
	* The size of the matrix A of the TwoBlockMatrix
	*/
	unsigned int
	N;

	/**
	 * UMFPACK symbolic decomposition of A
	 */
	void*
	symbolic_decomposition = nullptr;

	/**
	 * UMFPACK numeric decomposition of A
	 */
	void*
	numeric_decomposition = nullptr;

	/**
	 * Array saying which row starts where in BlockSolverWrapperUMFPACK2::Ai
	 */
	std::vector<SuiteSparse_long>
	Ap;

	/**
	 * Array storing the column indices of non-zero entries of A
	 */
	std::vector<SuiteSparse_long>
	Ai;

	/**
	 * Array storing the numerical values of the non-zero entries of the matrix A
	 */
	std::vector<double>
	Ax;

	/**
	 * Sets the data structures storing the matrix up (BlockSolverWrapperUMFPACK2::Ap, BlockSolverWrapperUMFPACK2::Ai, BlockSolverWrapperUMFPACK2::Ax).
	 *
	 * Also sets up the control parameters for the UMFPACK routines
	 *
	 * @param[in]	matrix		The matrix A
	 */
	void
	initialize_matrix(const SparseMatrix<double>& matrix);

	/**
	 * This function analyses the given matrix and sets up BlockSolverWrapperUMFPACK2::symbolic_decomposition.
	 * Whether this function does anything depends upon the value of BlockSolverWrapperUMFPACK2::analyze.
	 */
	void
	analyze_matrix();

	/**
	 * This function factorizes the matrix. The result is available from BlockSolverWrapperUMFPACK2::numeric_decomposition.
	 */
	void
	factorize_matrix();

	/**
	 * Multiply the inverse of A by @p f and store the result into @p x. I.e., solve Ax=f for x.
	 * This uses the factors computed by BlockSolverWrapperUMFPACK2::factorize_matrix().
	 *
	 * @param[out]	x	solution
	 * @param[in]	f	rhs
	 */
	void
	vmult(	Vector<double>& 		x,
			const Vector<double>&	f);

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
	unsigned int
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
			const bool											symmetric = false);

	DeclException2(	ExcUMFPACKError,
					std::string, int,
					<< "UMFPACK routine " << arg1 << " returned error status " << arg2 << ".");

};

#endif // DEAL_II_WITH_UMFPACK

// TODO: internally UMFPACK interface of deal.II is used to invert the matrix related to the second block, replace this by something more efficient!
#ifdef GALERKIN_TOOLS_WITH_MA57
#ifdef DEAL_II_WITH_UMFPACK

/** MA57 prototypes */
extern "C" void ma57id_(double*, int*);
extern "C" void ma57ad_(int*, int*, int*, int*, int*, int*, int*, int*, int*, double*);
extern "C" void ma57bd_(int*, int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, double*, int*, double*);
extern "C" void ma57cd_(int*, int*, double*, int*, int*, int*, int*, double*, int*, double*, int*, int*, int*, int*);
extern "C" void ma57dd_(int*, int*, int*, double*, int*, int*, double*, int*, int*, int*, double*, double*, double*, double*, int*, int*, double*, int*, double*);

/**
 * A SolverWrapper for the MA57 solver.
 *
 * @note MA57 has its own license, independent of that of GalerkinTools. If you
 * want to use MA57 you have to accept that license and obtain your own copy of MA57.
 *
 * Code partially copied over from deal.II
 */
class BlockSolverWrapperMA57 : public SolverWrapper<dealii::Vector<double>, dealii::BlockVector<double>, TwoBlockMatrix<dealii::SparseMatrix<double>>, TwoBlockSparsityPattern>
{
private:

	/**
	 * double control parameters of MA57
	 */
	std::vector<double>
	CNTL;

	/**
	 * integer control parameters of MA57
	 */
	std::vector<int>
	ICNTL;

	/**
	* The size of the matrix A of the TwoBlockMatrix
	*/
	int
	N;

	/**
	* The number of entries needed to define A (entries in upper triangular part of A)
	*/
	int
	NE;

	/**
	 * row indices of entries (attention: Fortran indexing->starts from 1)
	 */
	std::vector<int>
	IRN;

	/**
	 * column indices of entries (attention: Fortran indexing->starts from 1)
	 */
	std::vector<int>
	JCN;

	/**
	 * Array storing the numerical values of the matrix entrie (A(k) corresponds to row IRN(k), ICN(k))
	 */
	std::vector<double>
	A;

	/**
	 * Array for internal use of MA57
	 */
	std::vector<int>
	KEEP;

	/**
	 * length of KEEP
	 */
	int
	LKEEP;

	/**
	 * Workspace array for MA57
	 */
	std::vector<int>
	IWORK;

	/**
	 * another Workspace array for MA57
	 */
	std::vector<double>
	WORK;

	/**
	 * size if WORK
	 */
	int
	LWORK;

	/**
	 * Info array for MA57
	 */
	std::vector<int>
	INFO;

	/**
	 * RINFO array for MA57
	 */
	std::vector<double>
	RINFO;

	/**
	 * Factors of the matrix
	 */
	std::vector<double>
	FACT;

	/**
	 * entries in factors of the matrix
	 */
	int
	LFACT;

	/**
	 * Indexing information of factors of the matrix
	 */
	std::vector<int>
	IFACT;

	/**
	 * entries in IFACT
	 */
	int
	LIFACT;

	/**
	 * Sets the data structures storing the matrix up (BlockSolverWrapperMA57::Ap, BlockSolverWrapperMA57::Ai, BlockSolverWrapperMA57::Ax).
	 *
	 * Also sets up the control parameters for the MA57 routines
	 *
	 * @param[in]	matrix		The matrix A
	 */
	void
	initialize_matrix(	const SparseMatrix<double>& matrix);

	/**
	 * This function analyses the given matrix.
	 * Whether this function does anything depends upon the value of BlockSolverWrapperMA57::analyze.
	 */
	void
	analyze_matrix();

	/**
	 * This function factorizes the matrix.
	 */
	void
	factorize_matrix();

	/**
	 * Multiply the inverse of A by @p f and store the result into @p x. I.e., solve Ax=f for x.
	 * This uses the factors computed by BlockSolverWrapperMA57::factorize_matrix().
	 *
	 * @param[out]	x	solution
	 * @param[in]	f	rhs
	 */
	void
	vmult(	Vector<double>& 		x,
			const Vector<double>&	f);

public:

	/**
	 * Destructor.
	 */
	~BlockSolverWrapperMA57();

	/**
	 * Key indicating when to analyze the matrix structure:<br>
	 * 0 - before each factorization (default)<br>
	 * 1 - only before the next factorization (afterwards, BlockSolverWrapperMA57::analyze will be set to 2)<br>
	 * >=2 - do not recompute
	 *
	 * @warning		The user must ensure that BlockSolverWrapperMA57::analyze=0 or BlockSolverWrapperMA57::analyze=1
	 * 				whenever the sparsity pattern of the matrix has changed (or during the first call). Currently, no internal checking is performed.
	 *
	 * @warning		The scaling of the matrix is only recomputed if BlockSolverWrapperMA57::analyze <2. Otherwise, the most recently computed scaling will be used.
	 */
	unsigned int
	analyze = 0;

	/**
	 * diagnostic printing:<br>
	 * 0: no output<br>
	 * 1: print output
	 */
	unsigned int
	print_level = 0;

	/**
	 * ordering method (key corresponding to ICNTL(6) of MA57):<br>
	 * 0: AMD ordering using MC47 (without the detection of dense rows)<br>
	 * 2: AMD ordering using MC47<br>
	 * 3: Minimum degree ordering using MA27<br>
	 * 4: nested dissection (METIS)<br>
	 * 5: Automatic choice
	 * >=6: currently equivalent to 5
	 */
	unsigned int
	ordering_method = 5;

	/**
	 * indicates whether to use iterative refinement
	 */
	bool
	use_iterative_refinement = true;

	/**
	 * whether to ignore zeros in the matrix
	 */
	bool
	ignore_zeros = false;

	/**
	 * @copydoc SolverWrapper::solve
	 */
	virtual
	void
	solve(	const TwoBlockMatrix<dealii::SparseMatrix<double>>&	K_stretched,
			dealii::Vector<double>&								solution,
			const  dealii::BlockVector<double>&					f_stretched,
			const bool											symmetric = false);

	DeclException2(	ExcMA57Error,
					std::string, int,
					<< "MA57 routine " << arg1 << " returned error status " << arg2 << ".");

};

#endif // DEAL_II_WITH_UMFPACK
#endif // GALERKIN_TOOLS_WITH_MA57

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* INCREMENTALFE_SOLVERWRAPPER_H_ */
