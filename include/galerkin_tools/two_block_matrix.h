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

#ifndef GALERKIN_TOOLS_TWO_BLOCK_MATRIX_H_
#define GALERKIN_TOOLS_TWO_BLOCK_MATRIX_H_

#include <vector>

#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#ifdef DEAL_II_WITH_MPI
#include <mpi/mpi.h>
#endif

#include <galerkin_tools/config.h>
#include <galerkin_tools/two_block_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * This class implements a block matrix of the form
 *	\f{equation*}
 *  \begin{pmatrix}
 *  	\boldsymbol{A} & \boldsymbol{B} \\
 *  	\boldsymbol{C}^\top & \boldsymbol{D}
 *  \end{pmatrix}
 *  \f}
 *
 *  The main difference to the deal.II block matrices is that internally \f$\boldsymbol{C}\f$ is stored instead of the transpose of it.
 *  This is useful for parallel applications where finite element related dofs are coupled to "global" unknowns which are not related to any
 *  mesh.
 *
 *  @tparam		MatrixType		The internal type of the matrix blocks
 */
template<class MatrixType>
class TwoBlockMatrix
{

protected:

	/**
	 * The matrix A
	 */
	MatrixType
	A;

	/**
	 * The matrix B
	 */
	MatrixType
	B;

	/**
	 * The matrix C
	 */
	MatrixType
	C;

	/**
	 * The matrix D
	 */
	MatrixType
	D;

	/**
	 * The sparsity pattern of A
	 */
	SparsityPattern
	sp_A;

	/**
	 * The sparsity pattern of B
	 */
	SparsityPattern
	sp_B;

	/**
	 * The sparsity pattern of C
	 */
	SparsityPattern
	sp_C;

	/**
	 * The sparsity pattern of D
	 */
	SparsityPattern
	sp_D;

	/**
	 * The dimension of A
	 */
	unsigned int
	block_0_size = 0;

	/**
	 * The dimension of D
	 */
	unsigned int
	block_1_size = 0;

	/**
	 * The total dimension
	 */
	unsigned int
	total_dimension = 0;

public:

	using value_type = double;

	/**
	 * Default constructor
	 */
	TwoBlockMatrix() = default;

	/**
	 * Constructor for the sequential case
	 *
	 * @param[in]	sp		The sparsity pattern for the matrix
	 *
	 */
	TwoBlockMatrix(	const TwoBlockSparsityPattern&	sp);

	/**
	 * Virtual destructor to make the class polymorphic and allow dynamic casts.
	 */
	virtual
	~TwoBlockMatrix() = default;

	/**
	 * Reinitialize the object
	 *
	 * @param[in]	sp		The sparsity pattern for the matrix
	 */
	void
	reinit(	const TwoBlockSparsityPattern&	sp);

	/**
	 * The total number of rows
	 */
	unsigned int
	m()
	const;

	/**
	 * The total number of columns
	 */
	unsigned int
	n()
	const;

	/**
	 * %Function adding to an element of the matrix
	 *
	 * @param[in]	i		row
	 *
	 * @param[in]	j		column
	 *
	 * @param[in]	value	the number to be added
	 */
	void
	add(const unsigned int 	i,
		const unsigned int 	j,
		const double 		value);

	/**
	 * Add an array of values given by @p values in the given global matrix row at columns specified by @p col_indices in the sparse matrix.
	 *
	 * @param[in]	row						The row to add to
	 *
	 * @param[in]	n_cols					The number of values
	 *
	 * @param[in]	col_indices				The column indices
	 *
	 * @param[in]	values					The values
	 *
	 * @param[in]	elide_zero_values		This only exists for compatibility with deal.II and does not affect what's going on in this function
	 *
	 * @param[in]	col_indices_are_sorted	This only exists for compatibility with deal.II and does not affect what's going on in this function
	 *
	 */
	void
	add(const unsigned int	row,
		const unsigned int	n_cols,
		const unsigned int*	col_indices,
		const double*		values,
		const bool			elide_zero_values = false,
		const bool			col_indices_are_sorted = false);

	/**
	 * %Function setting an element of the matrix
	 *
	 * @param[in]	i		row
	 *
	 * @param[in]	j		column
	 *
	 * @param[in]	value	the number to be added
	 */
	void
	set(const unsigned int 	i,
		const unsigned int 	j,
		const double 		value);

	/**
	 * read access to element
	 *
	 * @param[in]	i		row
	 *
	 * @param[in]	j		column
	 *
	 */
	double
	operator()(	const unsigned int i,
				const unsigned int j)
	const;

	/**
	 * Return an element of the matrix (if the entry does not exist in the sparsity pattern, 0 is returned)
	 */
	double
	el(	const unsigned int i,
		const unsigned int j)
	const;

	/**
	 * %Function setting all elements of the matrix to @p value
	 *
	 * @param[in]	value	Only 0 is allowed for the value.
	 *
	 */
	TwoBlockMatrix<MatrixType>&
	operator=(const double value);

	/**
	 * %Function compressing this matrix
	 *
	 * @param[in]	operation	The operation (typically add or insert)
	 *
	 */
	void
	compress(const VectorOperation::values operation);

	/**
	 * @return		const reference to A
	 */
	const MatrixType&
	get_A()
	const;

	/**
	 * @return		const reference to B
	 */
	const MatrixType&
	get_B()
	const;

	/**
	 * @return		const reference to C
	 */
	const MatrixType&
	get_C()
	const;

	/**
	 * @return		const reference to D
	 */
	const MatrixType&
	get_D()
	const;

	/**
	 * @return		reference to A
	 */
	MatrixType&
	get_A();

	/**
	 * @return		reference to B
	 */
	MatrixType&
	get_B();

	/**
	 * @return		reference to C
	 */
	MatrixType&
	get_C();

	/**
	 * @return		reference to D
	 */
	MatrixType&
	get_D();

	/**
	 * @return	TwoBlockMatrix::block_0_size
	 */
	unsigned int
	get_block_0_size()
	const;

	/**
	 * @return	TwoBlockMatrix::block_1_size
	 */
	unsigned int
	get_block_1_size()
	const;

};

#ifdef DEAL_II_WITH_MPI
namespace parallel
{

/**
 * The parallel version of ::TwoBlockMatrix
 *
 *  @tparam		MatrixType		The internal type of the matrix blocks
 */
template<class MatrixType>
class TwoBlockMatrix : public dealii::GalerkinTools::TwoBlockMatrix<MatrixType>
{

public:

	TwoBlockMatrix() = default;

	/**
	 * Constructor for the parallel case
	 *
	 * @param[in]	sp							The sparsity pattern for the matrix
	 *
	 * @param[in]	locally_owned_indices		The locally owned rows of TwoBlockMatrix::A
	 *
	 * @param[in]	mpi_communicator			The MPI communicator to be used
	 *
	 */
	TwoBlockMatrix(	const TwoBlockSparsityPattern&	sp,
					const IndexSet&					locally_owned_indices,
					const MPI_Comm					mpi_communicator = MPI_COMM_WORLD);

	/**
	 * Reinitialize the object
	 *
	 * @param[in]	sp							The sparsity pattern for the matrix
	 *
	 * @param[in]	locally_owned_indices		The locally owned rows of TwoBlockMatrix::A
	 *
	 * @param[in]	mpi_communicator			The MPI communicator to be used
	 */
	void
	reinit(	const TwoBlockSparsityPattern&	sp,
			const IndexSet&					locally_owned_indices,
			const MPI_Comm					mpi_communicator = MPI_COMM_WORLD);

	/**
	 * %Function setting all elements of the matrix to @p value
	 *
	 * @param[in]	value	Only 0 is allowed for the value.
	 *
	 */
	dealii::GalerkinTools::parallel::TwoBlockMatrix<MatrixType>&
	operator=(const double value);

	/**
	 * @return	The mpi communicator
	 */
	const MPI_Comm&
	get_communicator()
	const;

};

}
#endif //DEAL_II_WITH_MPI

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKIN_TOOLS_TWO_BLOCK_MATRIX_H_ */
