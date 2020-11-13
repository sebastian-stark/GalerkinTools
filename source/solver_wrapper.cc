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
		direct_solver.vmult(solution, f_stretched.block(1));
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
		PETScWrappers::MPI::Vector solution_petsc(f_stretched.block(1).locally_owned_elements(), comm);
		PETScWrappers::SparseDirectMUMPS solver(cn, comm);
		solver.set_symmetric_mode(symmetric);
		solver.solve(D, solution_petsc, f_stretched.block(1));
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
{
	//mpi communicator
	MPI_Comm comm = K_stretched.get_communicator();

	//matrix sub blocks
	const auto& K = K_stretched.get_A();

	//size of bottom diagonal block of stretched system
	const unsigned int size_w = K_stretched.get_block_1_size();
	(void)size_w;

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

#ifdef GALERKIN_TOOLS_WITH_PARDISO
#ifdef DEAL_II_WITH_UMFPACK

void
BlockSolverWrapperPARDISO::initialize_matrix(	const SparseMatrix<double>& matrix)
{

	// if matrix is already initialized and it will be analyzed again, clear the entire memory because otherwise the memory will be leaking
	if( (analyze < 2) && initialized)
	{
		int phase = -1;
		int mtype = get_matrix_type();
		int error = 0;
		if(initialized)
			pardiso(pt, &maxfct, &mnum, &mtype, &phase, &N, nullptr, Ap.data(), Ai.data(), nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error,  dparm);
		AssertThrow(error == 0, ExcPARDISOError("pardiso", error));
	}

	msglvl = print_level;

	Assert(matrix.m() == matrix.n(), ExcNotQuadratic());

	N = matrix.m();
	const int N_nonzero = matrix_type < 2 ? matrix.n_nonzero_elements() : (matrix.n_nonzero_elements() - N) / 2 + N;
	Ap.resize(N + 1);
	Ai.resize(N_nonzero);
	Ax.resize(N_nonzero);

	Ap[0] = 1;
	unsigned int elements_in_row;
	if(matrix_type < 2)
	{
		for (int row = 0; row < N; ++row)
			Ap[row + 1] = Ap[row] + matrix.get_row_length(row);
		Assert(Ap.back() - 1 == N_nonzero, ExcMessage("Either this is a bug in the implementation of the PARDISO interface or the matrix supplied is not associated with a symmetric sparsity pattern!"));

		std::vector<int> row_pointers = Ap;
		for (int row = 0; row < N; ++row)
		{
			for(auto p = matrix.begin(row); p != matrix.end(row); ++p)
			{
				Ai[row_pointers[row] - 1] = p->column() + 1;
				Ax[row_pointers[row] - 1] = std::real(p->value());
				++row_pointers[row];
			}
			Assert(Ap[row + 1] == row_pointers[row], ExcMessage("Internal error - this indicates a bug!"));
		}
	}
	else
	{
		for(int row = 0; row < N; ++row)
		{
			elements_in_row = matrix.get_row_length(row);
			for(auto p = matrix.begin(row); p != matrix.end(row); ++p)
			{
				if(p->column() < p->row())
					elements_in_row--;
			}
			Ap[row + 1] = Ap[row] + elements_in_row;
		}
		Assert(Ap.back() - 1 == N_nonzero, ExcMessage("Either this is a bug in the implementation of the PARDISO interface or the matrix supplied is not associated with a symmetric sparsity pattern!"));

		std::vector<int> row_pointers = Ap;
		for (int row = 0; row < N; ++row)
		{
			for(auto p = matrix.begin(row); p != matrix.end(row); ++p)
			{
				if(p->column() >= p->row())
				{
					Ai[row_pointers[row] - 1] = p->column() + 1;
					Ax[row_pointers[row] - 1] = std::real(p->value());
					++row_pointers[row];
				}
			}
			Assert(Ap[row + 1] == row_pointers[row], ExcMessage("Internal error - this indicates a bug!"));
		}

	}

	// order entries (this swaps the diagonal element, which is stored first in deal.II, into the right position)
	for(int row = 0; row < N; ++row)
	{
		int cursor = Ap[row] - 1;
		while( ( cursor < Ap[row + 1] - 1 - 1 ) && ( Ai[cursor] > Ai[cursor + 1] ) )
		{
			std::swap(Ai[cursor], Ai[cursor + 1]);
			std::swap(Ax[cursor], Ax[cursor + 1]);
			++cursor;
		}
	}

	// set up control parameters
	if(analyze < 2)
	{
		int mtype = get_matrix_type();
		int solver = 0;
		int error = 0;
		initialized = true;
		pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);
		AssertThrow(error == 0, ExcPARDISOError("pardisoinit", error));

		if(!use_defaults)
		{
			Assert((ordering_method == 0) || (ordering_method == 2), ExcMessage("Only values 0 and 2 allowed for ordering method"));
			iparm[1] = ordering_method;

			Assert((apply_scaling == 0) || (apply_scaling == 1), ExcMessage("Only values 0 and 1 allowed for apply_scaling"));
			iparm[10] = iparm[12] = apply_scaling;
			if(matrix_type == 2)
			{
				Assert((pivoting_method == 0) || (pivoting_method == 1), ExcMessage("Only values 0 and 1 allowed for pivoting_method"));
				iparm[20] = pivoting_method;
			}
		    iparm[7] = n_iterative_refinements;
		    iparm[9] = 10;
		}

		// determine number of processors
		char* var = getenv("OMP_NUM_THREADS");
		int num_procs;
		if(var != NULL)
			sscanf(var, "%d", &num_procs);
		else
			Assert(false, ExcMessage("You have to export an environment variable OMP_NUM_THREADS in order to specify how many threads are allowed for PARDISO!"));
		iparm[2]  = num_procs;

/*		if(print_level == 1)
		{
			pardiso_printstats (&mtype, &N, Ax.data(), Ap.data(), Ai.data(), &nrhs, nullptr, &error);
			AssertThrow(error == 0, ExcPARDISOError("pardiso_printstats", error));
		}*/
	}

	return;
}

void
BlockSolverWrapperPARDISO::analyze_matrix()
{
	if(analyze < 2)
	{
		int phase = 11;
		int mtype = get_matrix_type();

		int error = 0;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &N, Ax.data(), Ap.data(), Ai.data(), nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error, dparm);
		AssertThrow(error == 0, ExcPARDISOError("pardiso", error));
	}

	// Matrix was only analyzed in this step
	if(analyze == 1)
		analyze = 2;

	return;
}


void
BlockSolverWrapperPARDISO::factorize_matrix()
{
	int phase = 22;
	int mtype = get_matrix_type();

	int error = 0;
	pardiso (pt, &maxfct, &mnum, &mtype, &phase, &N, Ax.data(), Ap.data(), Ai.data(), nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error,  dparm);
	if(error == -2)
		cout << "Ran out of memory, memory required=" << iparm[16]/1000000.0 << endl;
	AssertThrow(error == 0, ExcPARDISOError("pardiso", error));

	return;
}

void
BlockSolverWrapperPARDISO::vmult(	Vector<double>& 		x,
									const Vector<double>&	f )
{
	// make sure that some kind of factorize() call has happened before
	Assert(Ap.size() != 0, ExcNotInitialized());
	Assert(Ai.size() != 0, ExcNotInitialized());
	Assert(Ai.size() == Ax.size(), ExcNotInitialized());

	x.reinit(f.size());

	int phase = 33;
    int mtype = get_matrix_type();

    int error = 0;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &N, Ax.data(), Ap.data(), Ai.data(), nullptr, &nrhs, iparm, &msglvl, const_cast<double*>(f.data()), x.data(), &error,  dparm);

	/**
	 * check residual (TODO: it seems to be a bug in PARDISO that it does not complain about large residuals)
	 */
	double normb, normr;
	res.reinit(N, true);
	pardiso_residual (&mtype , &N, Ax.data(), Ap.data() , Ai.data(), const_cast<double*>(f.data()), x.data(), res.data(), &normb , &normr);
	if( normr >= res_max)
		cout << "The solution of the PARDISO solver is inaccurate! (Residual = " << normr << ")" << endl;
	AssertThrow(error == 0, ExcPARDISOError("pardiso", error));

	//AssertThrow(normr < res_max, ExcPARDISORes(normr));

	return;
}

int
BlockSolverWrapperPARDISO::get_matrix_type()
{
	switch(matrix_type)
	{
		case 0:
			return 11;
		case 1:
			return 1;
		case 2:
			return -2;
		default:
			return 2;
	}
}

BlockSolverWrapperPARDISO::~BlockSolverWrapperPARDISO()
{
	int phase = -1;
	int mtype = get_matrix_type();
	int error = 0;
	if(initialized)
		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &N, nullptr, Ap.data(), Ai.data(), nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error,  dparm);
	AssertThrow(error == 0, ExcPARDISOError("pardiso", error));
}

void
BlockSolverWrapperPARDISO::solve(const TwoBlockMatrix<SparseMatrix<double>>&	K_stretched,
		 	 	 	 	 		Vector<double>&									solution,
								const BlockVector<double>&						f_stretched,
								const bool 										/*symmetric*/)
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

	// initialize solver, analyze and factorize
	if(size_f > 0)
	{
		initialize_matrix(K);
		analyze_matrix();
		factorize_matrix();
	}

	// if one of the blocks is zero sized, solve directly with the other and return
	if(size_f == 0)
	{
		SparseDirectUMFPACK direct_solver_D;
		direct_solver_D.initialize(D);
		direct_solver_D.vmult(solution, f_stretched.block(0));
		return;
	}
	else if(size_w == 0)
	{
		vmult(solution, f_stretched.block(0));
		return;
	}

	//vector sub blocks
	const auto& f = f_stretched.block(0);
	const auto& w = f_stretched.block(1);

	//step 0 (u_0 = inv(K) * f)
	Vector<double> u_0(size_f);
	vmult(u_0, f);

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
		const auto L_U_n_ = L_U_n;
		vmult(L_U_n, L_U_n_);
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

	return;

}

#endif // GALERKIN_TOOLS_WITH_PARDISO
#endif // DEAL_II_WITH_UMFPACK

#ifdef DEAL_II_WITH_UMFPACK

void
BlockSolverWrapperUMFPACK2::initialize_matrix(const SparseMatrix<double>& matrix)
{
	// Code partially copied over from deal.II library

	Assert(matrix.m() == matrix.n(), ExcNotQuadratic());

	N = matrix.m();
	Ap.resize(N + 1);
	Ai.resize(matrix.n_nonzero_elements());
	Ax.resize(matrix.n_nonzero_elements());

	// copy data into the data structures (note that deal.II uses a slightly different format to store the matrix, which complicates matters a bit)
	Ap[0] = 0;
	for (unsigned int row = 1; row <= N; ++row)
		Ap[row] = Ap[row - 1] + matrix.get_row_length(row - 1);

	std::vector<SuiteSparse_long> row_pointers = Ap;
	for (unsigned int row = 0; row < matrix.m(); ++row)
	{
		for(auto p = matrix.begin(row); p != matrix.end(row); ++p)
		{
			Ai[row_pointers[row]] = p->column();
			Ax[row_pointers[row]] = std::real(p->value());
			++row_pointers[row];
		}
	}

	for(unsigned int row = 0; row < matrix.m(); ++row)
	{
		SuiteSparse_long cursor = Ap[row];
		while( ( cursor < Ap[row + 1] - 1 ) && ( Ai[cursor] > Ai[cursor + 1] ) )
		{
			std::swap(Ai[cursor], Ai[cursor + 1]);
			std::swap(Ax[cursor], Ax[cursor + 1]);
			++cursor;
		}
	}

	// set up control parameters
	umfpack_dl_defaults(control.data());

	return;
}

void
BlockSolverWrapperUMFPACK2::analyze_matrix()
{

	if(analyze < 2)
	{
		// throw away old symbolic decomposition if it exists
		if (symbolic_decomposition != nullptr)
		{
			umfpack_dl_free_symbolic(&symbolic_decomposition);
			symbolic_decomposition = nullptr;
		}

		// analyze the matrix

		const int status = umfpack_dl_symbolic(N, N, Ap.data(), Ai.data(), Ax.data(), &symbolic_decomposition, control.data(), info.data());
		AssertThrow(status == UMFPACK_OK, ExcUMFPACKError("umfpack_dl_symbolic", status));
	}

	// Matrix was only analyzed in this step
	if(analyze == 1)
		analyze = 2;

	return;
}

void
BlockSolverWrapperUMFPACK2::factorize_matrix()
{
	// throw away old numeric decomposition if it exists
	if (numeric_decomposition != nullptr)
	{
		umfpack_dl_free_numeric(&numeric_decomposition);
		numeric_decomposition = nullptr;
	}

	const int status = umfpack_dl_numeric(Ap.data(), Ai.data(), Ax.data(), symbolic_decomposition,  &numeric_decomposition, control.data(), info.data());
	AssertThrow(status == UMFPACK_OK, ExcUMFPACKError("umfpack_dl_numeric", status));

	// print report of factorization
	if(print_level > 0)
		cout << "*** UMFPACK REPORT OF FACTORIZATION ***" << endl << endl;
	control[UMFPACK_PRL] = print_level;
	umfpack_dl_report_info(control.data(), info.data());

	return;
}

void
BlockSolverWrapperUMFPACK2::vmult(	Vector<double>& 		x,
									const Vector<double>&	f )
{
	// make sure that some kind of factorize() call has happened before
	Assert(Ap.size() != 0, ExcNotInitialized());
	Assert(Ai.size() != 0, ExcNotInitialized());
	Assert(Ai.size() == Ax.size(), ExcNotInitialized());

	x.reinit(f.size());

	// solve the system. note that since UMFPACK wants compressed column
	// storage instead of the compressed row storage format we use in
	// deal.II's SparsityPattern classes, we solve for UMFPACK's A^T instead

	const int status = umfpack_dl_solve(UMFPACK_At, Ap.data(), Ai.data(), Ax.data(), x.begin(), f.begin(), numeric_decomposition, control.data(), info.data());
	AssertThrow(status == UMFPACK_OK, ExcUMFPACKError("umfpack_dl_solve", status));

	return;
}

BlockSolverWrapperUMFPACK2::~BlockSolverWrapperUMFPACK2()
{
	if (symbolic_decomposition != nullptr)
	{
		umfpack_dl_free_symbolic(&symbolic_decomposition);
		symbolic_decomposition = nullptr;
	}

	if (numeric_decomposition != nullptr)
	{
		umfpack_dl_free_numeric(&numeric_decomposition);
		numeric_decomposition = nullptr;
	}
}

void
BlockSolverWrapperUMFPACK2::solve(	const TwoBlockMatrix<SparseMatrix<double>>& K_stretched,
		 	 	 	 	 			Vector<double>&								solution,
									const BlockVector<double>&					f_stretched,
									const bool 									/*symmetric*/)
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

	// initialize solver, analyze and factorize
	if(size_f > 0)
	{
		initialize_matrix(K);
		analyze_matrix();

/*		FILE* printout = fopen ("K_2.dat","w");
		for (unsigned int row = 0; row < K.m(); ++row)
		{
			for(auto p = K.begin(row); p != K.end(row); ++p)
			{
				fprintf(printout, "%i %i %- 1.16e\n", p->row(), p->column(), std::real(p->value()));
			}
		}
		fclose(printout);
*/

		factorize_matrix();
	}

	// if one of the blocks is zero sized, solve directly with the other and return
	if(size_f == 0)
	{
		SparseDirectUMFPACK direct_solver_D;
		direct_solver_D.initialize(D);
		direct_solver_D.vmult(solution, f_stretched.block(0));
		return;
	}
	else if(size_w == 0)
	{
		vmult(solution, f_stretched.block(0));
		// throw away numeric decomposition
		if (numeric_decomposition != nullptr)
		{
			umfpack_dl_free_numeric(&numeric_decomposition);
			numeric_decomposition = nullptr;
		}

		return;
	}

	//vector sub blocks
	const auto& f = f_stretched.block(0);
	const auto& w = f_stretched.block(1);

	//step 0 (u_0 = inv(K) * f)
	Vector<double> u_0(size_f);
	vmult(u_0, f);

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
		const auto L_U_n_ = L_U_n;
		vmult(L_U_n, L_U_n_);
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

	// throw away numeric decomposition
	if (numeric_decomposition != nullptr)
	{
		umfpack_dl_free_numeric(&numeric_decomposition);
		numeric_decomposition = nullptr;
	}

	return;
}

#endif // DEAL_II_WITH_UMFPACK

#ifdef GALERKIN_TOOLS_WITH_MA57
#ifdef DEAL_II_WITH_UMFPACK

void
BlockSolverWrapperMA57::initialize_matrix(const SparseMatrix<double>& matrix)
{
	// Code partially copied over from deal.II library

	Assert(matrix.m() == matrix.n(), ExcNotQuadratic());

	N = matrix.m();
	NE = (matrix.n_nonzero_elements() - N) / 2 + N;

	IRN.resize(NE);
	JCN.resize(NE);
	A.resize(NE);

	int counter = 0;
	for (unsigned int row = 0; row < matrix.m(); ++row)
	{
		for(auto p = matrix.begin(row); p != matrix.end(row); ++p)
		{
			if(p->column() >= p->row())
			{
				if(!ignore_zeros)
				{
					IRN[counter] = p->row() + 1;
					JCN[counter] = p->column() + 1;
					A[counter] = std::real(p->value());
					++counter;
				}
				else
				{
					if(std::real(p->value()) != 0.0)
					{
						IRN[counter] = p->row() + 1;
						JCN[counter] = p->column() + 1;
						A[counter] = std::real(p->value());
						++counter;
					}
				}
			}
		}
	}
	if(ignore_zeros)
	{
		NE = counter;
		IRN.resize(NE);
		JCN.resize(NE);
		A.resize(NE);
	}

	CNTL.resize(5);
	ICNTL.resize(20);
	ma57id_(CNTL.data(), ICNTL.data());

	if(print_level == 1)
		ICNTL[4] = 3;
	Assert(ordering_method != 1, ExcMessage("ordering_method is not allowed to be 1!"));
	ICNTL[5] = ordering_method;

	return;
}

void
BlockSolverWrapperMA57::analyze_matrix()
{

	if(analyze < 2 || ignore_zeros)
	{
		INFO.resize(40);
		RINFO.resize(20);
		LKEEP = 5*N + NE + std::max(N,NE) + 42 + 2*N;
		KEEP.resize(LKEEP);
		IWORK.resize(5*N);
		LWORK = 5*N;
		WORK.resize(LWORK);

		ma57ad_(&N, &NE, IRN.data(), JCN.data(), &LKEEP, KEEP.data(), IWORK.data(), ICNTL.data(), INFO.data(), RINFO.data());
		AssertThrow(INFO[0] >= 0, ExcMA57Error("ma57ad_", INFO[0]));

		// allocate memory for factors; be conservative here -> allow for a lot of extra memory for pivoting for the indefinite case
		LFACT = INFO[8] * 2;
		LIFACT = INFO[9] * 2;
		FACT.resize(LFACT);
		IFACT.resize(LIFACT);
	}

	// Matrix was only analyzed in this step
	if(analyze == 1)
		analyze = 2;

	return;
}

void
BlockSolverWrapperMA57::factorize_matrix()
{

	ma57bd_(&N, &NE, A.data(), FACT.data(), &LFACT, IFACT.data(), &LIFACT, &LKEEP, KEEP.data(), IWORK.data(), ICNTL.data(), CNTL.data(), INFO.data(), RINFO.data());
	AssertThrow(INFO[0] >= 0, ExcMA57Error("ma57bd_", INFO[0]));

	return;
}

void
BlockSolverWrapperMA57::vmult(	Vector<double>& 		x,
								const Vector<double>&	f )
{

	int NRHS = 1;
	int JOB = 0;
	if(!use_iterative_refinement)
	{
		x = f;
		ma57cd_(&JOB, &N, FACT.data(), &LFACT, IFACT.data(), &LIFACT, &NRHS, x.data(), &N, WORK.data(), &LWORK, IWORK.data(), ICNTL.data(), INFO.data());

		AssertThrow(INFO[0] >= 0, ExcMA57Error("ma57cd_", INFO[0]));
	}
	else
	{
		x.reinit(f.size());
		Vector<double> RESID(f.size());
		ma57dd_(&JOB, &N, &NE, A.data(), IRN.data(), JCN.data(), FACT.data(), &LFACT, IFACT.data(), &LIFACT, const_cast<double*>(f.data()), x.data(), RESID.data(), WORK.data(), IWORK.data(), ICNTL.data(), CNTL.data(), INFO.data(), RINFO.data());
		AssertThrow(INFO[0] >= 0, ExcMA57Error("ma57dd_", INFO[0]));
	}

	AssertThrow(INFO[0] >= 0, ExcMA57Error("ma57cd_", INFO[0]));

	return;
}



void
BlockSolverWrapperMA57::solve(	const TwoBlockMatrix<SparseMatrix<double>>& K_stretched,
		 	 	 	 	 			Vector<double>&							solution,
									const BlockVector<double>&				f_stretched,
									const bool 								/*symmetric*/)
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

	// initialize solver, analyze and factorize
	if(size_f > 0)
	{
		initialize_matrix(K);
		analyze_matrix();
		factorize_matrix();
	}

	// if one of the blocks is zero sized, solve directly with the other and return
	if(size_f == 0)
	{
		SparseDirectUMFPACK direct_solver_D;
		direct_solver_D.initialize(D);
		direct_solver_D.vmult(solution, f_stretched.block(0));
		return;
	}
	else if(size_w == 0)
	{
		vmult(solution, f_stretched.block(0));
		return;
	}

	//vector sub blocks
	const auto& f = f_stretched.block(0);
	const auto& w = f_stretched.block(1);

	//step 0 (u_0 = inv(K) * f)
	Vector<double> u_0(size_f);
	vmult(u_0, f);

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
		const auto L_U_n_ = L_U_n;
		vmult(L_U_n, L_U_n_);
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

	return;
}

BlockSolverWrapperMA57::~BlockSolverWrapperMA57()
{
}

#endif // DEAL_II_WITH_UMFPACK
#endif // GALERKIN_TOOLS_WITH_MA57


template class SolverWrapper<Vector<double>, Vector<double>, SparseMatrix<double>, SparsityPattern>;
template class SolverWrapper<Vector<double>, BlockVector<double>, TwoBlockMatrix<SparseMatrix<double>>, TwoBlockSparsityPattern>;
#ifdef DEAL_II_WITH_PETSC
#ifdef DEAL_II_WITH_MPI
	template class SolverWrapper<LinearAlgebra::distributed::Vector<double>, PETScWrappers::MPI::BlockVector, dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>, TwoBlockSparsityPattern>;
#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_WITH_MPI

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
