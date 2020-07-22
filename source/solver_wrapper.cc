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

#include <galerkin_tools/solver_wrapper.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/petsc_precondition.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<class SolutionVectorType, class RHSVectorType, class MatrixType, class SparsityPatternType>
SolverWrapper<SolutionVectorType, RHSVectorType, MatrixType, SparsityPatternType>::~SolverWrapper()
{
	Assert(n_subscriptions() == 0, ExcMessage("You are about to destroy a SolverWrapper, which is currently in use! Make sure that all SolverWrapper objects live as least as long as the objects using them!"));
}

void
SolverWrapperUMFPACK::solve(const SparseMatrix<double>& K_stretched,
		 	 	 	 	 	Vector<double>&				solution,
							const Vector<double>&		f_stretched,
							const bool 					/*symmetric*/)
const
{
	SparseDirectUMFPACK direct_solver;
	direct_solver.initialize(K_stretched);
	direct_solver.vmult (solution, f_stretched);
}

void
BlockSolverWrapperUMFPACK::solve(	const TwoBlockMatrix<SparseMatrix<double>>& K_stretched,
		 	 	 	 	 			Vector<double>&				solution,
									const BlockVector<double>&	f_stretched,
									const bool 					/*symmetric*/)
const
{

	//matrix sub blocks
	const auto& K = K_stretched.get_A();
	const auto& L_U = K_stretched.get_B();
	const auto& L_V = K_stretched.get_C();
	const auto& D = K_stretched.get_D();

	//size of top diagonal block of stretched system
	const unsigned int size_f = K_stretched.get_block_0_size();
	//size of bottom diagonal block of stretched system
	const unsigned int size_w = K_stretched.get_block_1_size();

	//initialize solver
	SparseDirectUMFPACK direct_solver;
	if(size_f == 0)
	{
		direct_solver.initialize(D);
		direct_solver.vmult(solution, f_stretched.block(0));
		return;
	}
	else if(size_w == 0)
	{
		direct_solver.initialize(K);
		direct_solver.vmult(solution, f_stretched.block(0));
		return;
	}
	else
	{
		direct_solver.initialize(K);
	}

	//vector sub blocks
	const auto& f = f_stretched.block(0);
	const auto& w = f_stretched.block(1);

	//step 0 (u_0 = inv(K) * f)
	Vector<double> u_0(size_f);
	direct_solver.vmult(u_0, f);

	//step 1 (compute C = inv(K)*L_U)
	DynamicSparsityPattern dsp_C(size_f, size_w);
	for(unsigned int m = 0; m < size_f; ++m)
		for(unsigned int n = 0; n < size_w; ++n)
			dsp_C.add(m, n);
	SparsityPattern sp_C;
	sp_C.copy_from(dsp_C);
	SparseMatrix<double> C(sp_C);
	Vector<double> L_U_n(size_f);
	for(unsigned int n = 0; n < size_w; ++n)
	{
		for(unsigned int m = 0; m < size_f; ++m)
			L_U_n[m] = L_U.el(m, n);
		direct_solver.solve(L_U_n);
		for(unsigned int m = 0; m < size_f; ++m)
			C.set(m, n, L_U_n[m]);
	}

	//step 2 (compute F = L_V^T*C - D)
	DynamicSparsityPattern dsp_F(size_w, size_w);
	for(unsigned int m = 0; m < size_w; ++m)
		for(unsigned int n = 0; n < size_w; ++n)
			dsp_F.add(m, n);
	SparsityPattern sp_F;
	sp_F.copy_from(dsp_F);
	SparseMatrix<double> F(sp_F);
	L_V.Tmmult(F, C, Vector<double>(), false);
	for(unsigned int m = 0; m < size_w; ++m)
		for(unsigned int n = 0; n < size_w; ++n)
			F.add(m, n, -D.el(m, n));

	//step 3 (compute -lambda = inv(F) * (w - L_V^T*u_0))
	Vector<double> minus_lambda(size_w);
	L_V.Tvmult(minus_lambda, u_0);
	for(unsigned int m = 0; m < size_w; ++m)
		minus_lambda[m] = w[m] - minus_lambda[m];
	SparseDirectUMFPACK inv_F;
	inv_F.initialize(F);
	inv_F.solve(minus_lambda);

	//step 4 (compute u = u_0 - C*lambda)
	C.vmult_add(u_0, minus_lambda);

	//step 5 (transfer solution)
	for(unsigned int m = 0; m < size_f; ++m)
		solution[m] = u_0[m];
	for(unsigned int m = 0; m < size_w; ++m)
		solution[m + size_f] = -minus_lambda[m];
}

#ifdef DEAL_II_WITH_PETSC
#ifdef DEAL_II_PETSC_WITH_MUMPS
#ifdef DEAL_II_WITH_MPI

void
SolverWrapperPETSc::solve(	const parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>&	K_stretched,
							LinearAlgebra::distributed::Vector<double>&							solution,
							const PETScWrappers::MPI::BlockVector&								f_stretched,
							const bool															symmetric)
const
{

	//mpi communicator
	MPI_Comm comm = K_stretched.get_communicator();

	//matrix sub blocks
	const auto& K = K_stretched.get_A();
	const auto& L_U = K_stretched.get_B();
	const auto& L_V = K_stretched.get_C();
	const auto& D = K_stretched.get_D();

	//size of top diagonal block of stretched system
	const unsigned int size_f = K_stretched.get_block_0_size();
	//size of bottom diagonal block of stretched system
	const unsigned int size_w = K_stretched.get_block_1_size();

	//solver control for MUMPS
	SolverControl cn;

	if(size_f == 0)
	{
		PETScWrappers::MPI::Vector solution_petsc(f_stretched.block(0).locally_owned_elements(), comm);
		PETScWrappers::SparseDirectMUMPS solver(cn, comm);
		solver.set_symmetric_mode(symmetric);
		solver.solve(D, solution_petsc, f_stretched.block(0));
		for(const auto m : solution_petsc.locally_owned_elements())
			solution(m) = solution_petsc(m);
		solution.compress(VectorOperation::insert);
		return;
	}
	else if(size_w == 0)
	{
		PETScWrappers::MPI::Vector solution_petsc(f_stretched.block(0).locally_owned_elements(), comm);
		PETScWrappers::SparseDirectMUMPS solver(cn, comm);
		solver.set_symmetric_mode(symmetric);
		solver.solve(K, solution_petsc, f_stretched.block(0));
		for(const auto m : solution_petsc.locally_owned_elements())
			solution(m) = solution_petsc(m);
		solution.compress(VectorOperation::insert);
		return;
	}

	//vector sub blocks
	const auto& f = f_stretched.block(0);
	const auto& w = f_stretched.block(1);

	//locally owned indices of blocks
	const auto indices_block_0 = f.locally_owned_elements();
	const auto indices_block_1 = w.locally_owned_elements();

	//all indices of block 1
	IndexSet all_indices_block_1(size_w);
	all_indices_block_1.add_range(0, size_w);

	//step 0 (u_0 = inv(K) * f)
	PETScWrappers::MPI::Vector u_0(indices_block_0, comm);
	PETScWrappers::SparseDirectMUMPS solver(cn, comm);
	solver.set_symmetric_mode(symmetric);
	solver.solve(K, u_0, f);

	//step 1 (compute C = inv(K)*L_U)
	DynamicSparsityPattern dsp_C(size_f, size_w, indices_block_0);
	for(const auto& row : indices_block_0)
		for(const auto& col : all_indices_block_1)
			dsp_C.add(row, col);
	PETScWrappers::MPI::SparseMatrix C;
	C.reinit(indices_block_0, indices_block_1, dsp_C, comm);
	PETScWrappers::MPI::Vector L_U_n(indices_block_0, comm);
	PETScWrappers::MPI::Vector C_n(indices_block_0, comm);
	for(unsigned int n = 0; n < size_w; ++n)
	{
		for(const auto& m : indices_block_0)
			L_U_n(m) = L_U(m, n);
		L_U_n.compress(VectorOperation::insert);
		solver.solve(K, C_n, L_U_n);
		for(const auto& m : indices_block_0)
			C.set(m, n, C_n(m));
		C.compress(VectorOperation::insert);
	}

	//step 2 (compute F = L_V^T*C - D)
	DynamicSparsityPattern dsp_F(size_w, size_w, indices_block_1);
	for(const auto& m : indices_block_1)
		for(const auto& n : all_indices_block_1)
			dsp_F.add(m, n);
	PETScWrappers::MPI::SparseMatrix L_V_T_C, F;
	L_V_T_C.reinit(indices_block_1, indices_block_1, dsp_F, comm);
	F.reinit(indices_block_1, indices_block_1, dsp_F, comm);
	//todo: this seems to be a bug in PETSc -> look into it
	Assert(indices_block_0.n_elements() != 0, ExcMessage("SolverWrapperPETSc does currently not work if there is a processor not owning any elements!"));
	L_V.Tmmult(L_V_T_C, C);
	for(const auto& m : indices_block_1)
		for(const auto& n : all_indices_block_1)
			F.set(m, n, L_V_T_C.el(m, n) - D.el(m,n));
	F.compress(VectorOperation::insert);

	//step 3 (compute -lambda = inv(F) * (w - L_V^T*u_0))
	PETScWrappers::MPI::Vector w_minus_L_V_T_u_0(indices_block_1, comm);
	L_V.Tvmult(w_minus_L_V_T_u_0, u_0);
	for(const auto& m : indices_block_1)
		w_minus_L_V_T_u_0[m] = w[m] - w_minus_L_V_T_u_0[m];
	w_minus_L_V_T_u_0.compress(VectorOperation::insert);

	PETScWrappers::MPI::Vector minus_lambda(indices_block_1, comm);
	PETScWrappers::SparseDirectMUMPS solver_2(cn, comm);
	solver_2.set_symmetric_mode(symmetric);
	solver_2.solve(F, minus_lambda, w_minus_L_V_T_u_0);


	//step 4 (compute u = u_0 - C*lambda)
	C.vmult_add(u_0, minus_lambda);

	//step 5 (transfer solution)
	solution.zero_out_ghosts();
	for(const auto& m : indices_block_0)
		solution[m] = u_0[m];
	for(const auto& m : indices_block_1)
		solution[m + size_f] = -minus_lambda[m];
	solution.compress(VectorOperation::insert);
	solution.update_ghost_values();
}

#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_PETSC_WITH_MUMPS
#endif // DEAL_II_WITH_MPI


#ifdef DEAL_II_WITH_PETSC
#ifdef DEAL_II_WITH_MPI

void
SolverWrapperPETScIterative::solve(	const parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>&	K_stretched,
									LinearAlgebra::distributed::Vector<double>&							solution,
									const PETScWrappers::MPI::BlockVector&								f_stretched,
									const bool															symmetric)
const
{
	//mpi communicator
	MPI_Comm comm = K_stretched.get_communicator();

	//matrix sub blocks
	const auto& K = K_stretched.get_A();

	//size of bottom diagonal block of stretched system
	const unsigned int size_w = K_stretched.get_block_1_size();

	//solver control for MUMPS
	SolverControl cn;

	Assert(size_w == 0, ExcMessage("This is currently only implemented for systems without a second block!"));

	PETScWrappers::MPI::Vector solution_petsc(f_stretched.block(0).locally_owned_elements(), comm);

	SolverControl solver_control(solution.size(), 1e-12);
	PETScWrappers::SolverGMRES solver(solver_control, comm);

	PETScWrappers::PreconditionBoomerAMG preconditioner;
	PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
	data.symmetric_operator = symmetric;
    preconditioner.initialize(K, data);
    solver.solve(K,	solution_petsc, f_stretched.block(0), preconditioner);

    for(const auto m : solution_petsc.locally_owned_elements())
		solution(m) = solution_petsc(m);
	solution.compress(VectorOperation::insert);

	return;
}

#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_WITH_MPI



template class SolverWrapper<Vector<double>, Vector<double>, SparseMatrix<double>, SparsityPattern>;
template class SolverWrapper<Vector<double>, BlockVector<double>, TwoBlockMatrix<SparseMatrix<double>>, TwoBlockSparsityPattern>;
#ifdef DEAL_II_WITH_PETSC
#ifdef DEAL_II_WITH_MPI
	template class SolverWrapper<LinearAlgebra::distributed::Vector<double>, PETScWrappers::MPI::BlockVector, dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>, TwoBlockSparsityPattern>;
#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_WITH_MPI

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
