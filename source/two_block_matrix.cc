// --------------------------------------------------------------------------
// Copyright (C) 2019 by Sebastian Stark
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

#include <galerkin_tools/two_block_matrix.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/affine_constraints.templates.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<class MatrixType>
TwoBlockMatrix<MatrixType>::TwoBlockMatrix(	const TwoBlockSparsityPattern&	sp)
{
	reinit(sp);
}

template<class MatrixType>
void
TwoBlockMatrix<MatrixType>::reinit(	const TwoBlockSparsityPattern&	sp)
{
	block_0_size = sp.get_sp_A().n_rows();
	block_1_size = sp.get_sp_D().n_rows();
	total_dimension = block_0_size + block_1_size;

	auto A_ptr = dynamic_cast<SparseMatrix<double>*>(&A);
	auto B_ptr = dynamic_cast<SparseMatrix<double>*>(&B);
	auto C_ptr = dynamic_cast<SparseMatrix<double>*>(&C);
	auto D_ptr = dynamic_cast<SparseMatrix<double>*>(&D);

	Assert(A_ptr != nullptr, ExcMessage("The sequential version of TwoBlockMatrix is currently only implemented for SparseMatrix<double>!"));

	A_ptr->reinit(sp.get_sp_A());
	B_ptr->reinit(sp.get_sp_B());
	C_ptr->reinit(sp.get_sp_C());
	D_ptr->reinit(sp.get_sp_D());
}

template<class MatrixType>
unsigned int
TwoBlockMatrix<MatrixType>::m()
const
{
	return total_dimension;
}

template<class MatrixType>
unsigned int
TwoBlockMatrix<MatrixType>::n()
const
{
	return total_dimension;
}

template<class MatrixType>
void
TwoBlockMatrix<MatrixType>::add(const unsigned int 	i,
								const unsigned int 	j,
								const double 		value)
{
	if(i < block_0_size)
	{
		if(j < block_0_size)
			A.add(i, j, value);
		else
			B.add(i, j - block_0_size, value);
	}
	else
	{
		if(j < block_0_size)
			C.add(j, i - block_0_size, value);
		else
			D.add(i - block_0_size, j - block_0_size, value);
	}
}

template<class MatrixType>
void
TwoBlockMatrix<MatrixType>::add(const unsigned int	row,
								const unsigned int	n_cols,
								const unsigned int*	col_indices,
								const double*		values,
								const bool			/*elide_zero_values*/,
								const bool			/*col_indices_are_sorted*/)
{
	static std::vector<unsigned int> col_indices_0, col_indices_1;
	col_indices_0.resize(n_cols);
	col_indices_1.resize(n_cols);

	static std::vector<double> values_0, values_1;
	values_0.resize(n_cols);
	values_1.resize(n_cols);

	unsigned int elem_count_0 = 0;
	unsigned int elem_count_1 = 0;
	for(unsigned int col = 0; col < n_cols; ++col)
	{
		if(col_indices[col] < block_0_size)
		{
			col_indices_0[elem_count_0] = col_indices[col];
			values_0[elem_count_0] = values[col];
			++elem_count_0;
		}
		else
		{
			col_indices_1[elem_count_1] = col_indices[col] - block_0_size;
			values_1[elem_count_1] = values[col];
			++elem_count_1;
		}
	}
	col_indices_0.resize(elem_count_0);
	col_indices_1.resize(elem_count_1);
	values_0.resize(elem_count_0);
	values_1.resize(elem_count_1);

	if( row < block_0_size )

	{
		if(elem_count_0 > 0)
			A.add(row, elem_count_0, col_indices_0.data(), values_0.data());
		if(elem_count_1 > 0)
			B.add(row, elem_count_1, col_indices_1.data(), values_1.data());
	}
	else
	{
		for(unsigned int col = 0; col < elem_count_0; ++col)
			C.add(col_indices_0[col], row - block_0_size, values_0[col]);
		if(elem_count_1 > 0)
			D.add(row - block_0_size, elem_count_1, col_indices_1.data(), values_1.data());
	}
}

template<class MatrixType>
void
TwoBlockMatrix<MatrixType>::set(const unsigned int 	i,
								const unsigned int 	j,
								const double 		value)
{
	if(i < block_0_size)
	{
		if(j < block_0_size)
			A.set(i, j, value);
		else
			B.set(i, j - block_0_size, value);
	}
	else
	{
		if(j < block_0_size)
			C.set(j, i - block_0_size, value);
		else
			D.set(i - block_0_size, j - block_0_size, value);
	}
}

template<class MatrixType>
double
TwoBlockMatrix<MatrixType>::operator()(	const unsigned int i,
										const unsigned int j)
const
{
	if(i < block_0_size)
	{
		if(j < block_0_size)
			return A(i, j);
		else
			return B(i, j - block_0_size);
	}
	else
	{
		if(j < block_0_size)
			return C(j, i - block_0_size);
		else
			return D(i - block_0_size, j - block_0_size);
	}
}

template<class MatrixType>
double
TwoBlockMatrix<MatrixType>::el(	const unsigned int i,
								const unsigned int j)
const
{
	if(i < block_0_size)
	{
		if(j < block_0_size)
			return A.el(i, j);
		else
			return B.el(i, j - block_0_size);
	}
	else
	{
		if(j < block_0_size)
			return C.el(j, i - block_0_size);
		else
			return D.el(i - block_0_size, j - block_0_size);
	}
}

template<class MatrixType>
TwoBlockMatrix<MatrixType>&
TwoBlockMatrix<MatrixType>::operator=(const double value)
{
	(void)value;
	Assert(value == 0.0, ExcMessage("This operation is only allowed with value == 0.0"));
	if( (A.m() > 0) && (A.m() > 0))
		A = 0.0;
	if( (B.m() > 0) && (B.m() > 0))
		B = 0.0;
	if( (C.m() > 0) && (C.m() > 0))
		C = 0.0;
	if( (D.m() > 0) && (D.m() > 0))
		D = 0.0;
	return *this;
}

template<class MatrixType>
void
TwoBlockMatrix<MatrixType>::compress(const VectorOperation::values operation)
{
	A.compress(operation);
	B.compress(operation);
	C.compress(operation);
	D.compress(operation);
}

template<class MatrixType>
const MatrixType&
TwoBlockMatrix<MatrixType>::get_A()
const
{
	return A;
}

template<class MatrixType>
const MatrixType&
TwoBlockMatrix<MatrixType>::get_B()
const
{
	return B;
}

template<class MatrixType>
const MatrixType&
TwoBlockMatrix<MatrixType>::get_C()
const
{
	return C;
}

template<class MatrixType>
const MatrixType&
TwoBlockMatrix<MatrixType>::get_D()
const
{
	return D;
}

template<class MatrixType>
MatrixType&
TwoBlockMatrix<MatrixType>::get_A()

{
	return A;
}

template<class MatrixType>
MatrixType&
TwoBlockMatrix<MatrixType>::get_B()

{
	return B;
}

template<class MatrixType>
MatrixType&
TwoBlockMatrix<MatrixType>::get_C()
{
	return C;
}

template<class MatrixType>
MatrixType&
TwoBlockMatrix<MatrixType>::get_D()
{
	return D;
}

template<class MatrixType>
unsigned int
TwoBlockMatrix<MatrixType>::get_block_0_size()
const
{
	return block_0_size;
}

template<class MatrixType>
unsigned int
TwoBlockMatrix<MatrixType>::get_block_1_size()
const
{
	return block_1_size;
}

namespace parallel
{

#ifdef DEAL_II_WITH_MPI

template<class MatrixType>
TwoBlockMatrix<MatrixType>::TwoBlockMatrix(	const TwoBlockSparsityPattern&	dsp,
											const IndexSet&					locally_owned_indices,
											const MPI_Comm					mpi_communicator)
{
	reinit(dsp, locally_owned_indices, mpi_communicator);
}

template<class MatrixType>
void
TwoBlockMatrix<MatrixType>::reinit(	const TwoBlockSparsityPattern&	sp,
									const IndexSet&					locally_owned_indices,
									const MPI_Comm					mpi_communicator)
{
	this->block_0_size = sp.get_sp_A().n_rows();
	this->block_1_size = sp.get_sp_D().n_rows();
	this->total_dimension = this->block_0_size + this->block_1_size;

	const IndexSet indices_block_0 = locally_owned_indices.get_view(0, this->block_0_size);
	const IndexSet indices_block_1 = locally_owned_indices.get_view(this->block_0_size, this->total_dimension);

	if(this->block_0_size > 0)
		this->A.reinit(indices_block_0, indices_block_0, sp.get_sp_A(), mpi_communicator);
	if( (this->block_0_size > 0) && (this->block_1_size > 0))
	{
		this->B.reinit(indices_block_0, indices_block_1, sp.get_sp_B(), mpi_communicator);
		this->C.reinit(indices_block_0, indices_block_1, sp.get_sp_C(), mpi_communicator);
	}
	if(this->block_1_size > 0)
		this->D.reinit(indices_block_1, indices_block_1, sp.get_sp_D(), mpi_communicator);
}

template<class MatrixType>
dealii::GalerkinTools::parallel::TwoBlockMatrix<MatrixType>&
TwoBlockMatrix<MatrixType>::operator=(const double value)
{
	(void)value;
	Assert(value == 0.0, ExcMessage("This operation is only allowed with value == 0.0"));
	if( (this->A.m() > 0) && (this->A.m() > 0))
		this->A = 0.0;
	if( (this->B.m() > 0) && (this->B.m() > 0))
		this->B = 0.0;
	if( (this->C.m() > 0) && (this->C.m() > 0))
		this->C = 0.0;
	if( (this->D.m() > 0) && (this->D.m() > 0))
		this->D = 0.0;
	return *this;
}

template<class MatrixType>
const MPI_Comm&
TwoBlockMatrix<MatrixType>::get_communicator()
const
{
	return this->A.get_mpi_communicator();
}

#endif //DEAL_II_WITH_MPI

}


#ifdef DEAL_II_WITH_PETSC
	template class dealii::GalerkinTools::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>;
	template class dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>;
#endif //DEAL_II_WITH_PETSC

template class TwoBlockMatrix<SparseMatrix<double>>;

DEAL_II_NAMESPACE_CLOSE
GALERKIN_TOOLS_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

template
void
AffineConstraints<double>::distribute_local_to_global<GalerkinTools::TwoBlockMatrix<SparseMatrix<double>>>(	const FullMatrix<double>&,
																											const vector<AffineConstraints<double>::size_type>&,
																											const vector<AffineConstraints<double>::size_type>&,
																											GalerkinTools::TwoBlockMatrix<SparseMatrix<double>> &)
const;

template
void
AffineConstraints<double>::distribute_local_to_global<GalerkinTools::TwoBlockMatrix<SparseMatrix<double>>, BlockVector<double>>(const FullMatrix<double>&,
																																const Vector<double>&,
																																const vector<AffineConstraints<double>::size_type>&,
																																GalerkinTools::TwoBlockMatrix<SparseMatrix<double>>&,
																																BlockVector<double>&,
																																bool,
																																integral_constant<bool, false>)
const;

template
void
AffineConstraints<double>::distribute_local_to_global<dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>>(	const FullMatrix<double>&,
																																			const vector<AffineConstraints<double>::size_type>&,
																																			const vector<AffineConstraints<double>::size_type>&,
																																			dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>&)
const;

template
void
AffineConstraints<double>::distribute_local_to_global<dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>, dealii::PETScWrappers::MPI::BlockVector>(	const FullMatrix<double>&,
																																													const Vector<double>&,
																																													const vector<AffineConstraints<double>::size_type>&,
																																													dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>&,
																																													dealii::PETScWrappers::MPI::BlockVector&,
																																													bool,
																																													integral_constant<bool, false>)
const;

DEAL_II_NAMESPACE_CLOSE
