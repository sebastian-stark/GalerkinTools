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

#include <galerkin_tools/ldr.h>

#include <math.h>

#include <deal.II/base/exceptions.h>

#include <lapacke.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

int
Auxiliary::compute_ldr(	FullMatrix<double>& 			C,
						Vector<double>& 				D,
						vector<Vector<double>>& 		L,
						std::vector<Vector<double>>&	R)
{

	//setup
	int matrix_layout=LAPACK_ROW_MAJOR;
	char jobu='A';
	char jobvt='A';
	lapack_int m=C.size()[0];
	lapack_int n=C.size()[1];
	unique_ptr<double[]> a(new double[m*n]);
	int lda=m;
	unique_ptr<double[]> s(new double[m]);
	unique_ptr<double[]> u(new double[m*m]);
	lapack_int ldu=m;
	unique_ptr<double[]> vt(new double[m*m]);
	lapack_int ldvt=m;
	unique_ptr<double[]> superb(new double[m-1]);

	for(int i1=0; i1<m; i1++)
		for(int i2=0; i2<m; i2++)
			a[i1*m+i2]=C(i1,i2);

	Assert(n==m, ExcMessage("Matrix not square!"));

	//perform decomposition
	int error_code=LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a.get(), lda, s.get(), u.get(), ldu, vt.get(), ldvt, superb.get());
	Assert(error_code==0, ExcMessage("Something has went wrong in singular value decomposition LAPACKE_dgesvd!"));

	//transfer results into relevant data structures
	D.reinit(m);
	L.resize(m);
	R.resize(m);
	for(int i1=0; i1<m; i1++){
		L[i1].reinit(m);
		R[i1].reinit(m);
		D[i1]=s[i1];
		for(int i2=0; i2<m; i2++){
			L[i1][i2]=u[i1*m+i2];
			R[i1][i2]=vt[i2*m+i1];
		}
	}

	//scale entries of D to absolute values of 1 and ensure that L=R if C is symmetric
	for(int i2=0; i2<m; i2++){
		double sqrt_d=sqrt(D[i2]);
		D[i2]=1.0;
		for(int i1=0; i1<m; i1++){
			R[i1][i2]*=sqrt_d;
			L[i1][i2]*=sqrt_d;
		}
		double sp=0.0;
		for(int i1=0; i1<m; i1++)
			sp+=L[i1][i2]*R[i1][i2];
		if(sp<0.0){
			for(int i1=0; i1<m; i1++)
				R[i1][i2]*=-1.0;
			D[i2]=-1.0;
		}
	}

	return error_code;
}

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
