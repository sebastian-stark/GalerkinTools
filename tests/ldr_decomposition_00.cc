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

#include <iostream>
#include <vector>
#include <math.h>

#include <galerkin_tools/ldr.h>

using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;

void check(const bool sym_mode)
{
	unsigned int dim = 10;			//dimension

	FullMatrix<double> C(dim);
	Vector<double> D;
	vector<Vector<double>> L;
	vector<Vector<double>> R;

	//generate some matrix
	const unsigned int size = C.size()[0];
	for(unsigned int m = 0; m < size; ++m)
		for(unsigned int n = 0; n < size; ++n)
			C(m, n) = sin(sqrt(fabs(2*numbers::PI * double(m - n)/double(size)))) + sin(sqrt(fabs(2*numbers::PI * double(m + n)/double(size))));


	//symmetrize if required
	if(sym_mode)
		C.symmetrize();

	//perform decomposition
	Auxiliary::compute_ldr(C, D, L, R);

	//analyze result and check for errors
		FullMatrix<double> L_(dim);
	FullMatrix<double> D_(dim);
	FullMatrix<double> RT_(dim);
	FullMatrix<double> R_(dim);
	FullMatrix<double> D_RT(dim);
	FullMatrix<double> L_D_RT(dim);
	for(unsigned int m = 0; m < dim; ++m)
	{
		D_(m, m) = D(m);
		for(unsigned int n = 0; n < dim; ++n)
		{
			L_(m, n) = L[m][n];
			RT_(m, n) = R[n][m];
			R_(n, m) = RT_(m, n);
		}
	}
	D_.mmult(D_RT, RT_);
	L_.mmult(L_D_RT, D_RT);

	FullMatrix<double> L_D_RT_C(dim);
	for(unsigned int m = 0; m < dim; ++m)
		for(unsigned int n = 0; n < dim; ++n)
			L_D_RT_C(m, n) = L_D_RT(m, n)-C(m, n);

	cout << "L*D*R'-C" << endl;
	for(unsigned int m = 0; m < dim; ++m)
	{
		for(unsigned int n = 0; n < dim; ++n)
		{
			if(fabs(L_D_RT_C(m,n)) < 1e-14)
				L_D_RT_C(m, n) = 0.0;
			printf("% 10.8f ", L_D_RT_C(m, n));
		}
		cout << endl;
	}
	cout << endl;

	FullMatrix<double> U_V(dim);
	for(unsigned int m = 0; m < dim; ++m)
		for(unsigned int n=0; n<dim; ++n)
			U_V(m, n) = L_(m, n) - R_(m, n);

	cout << "L-R" << endl;
	for(unsigned int m = 0; m < dim; ++m)
	{
		for(unsigned int n = 0; n < dim; ++n)
		{
			if(fabs(U_V(m, n)) < 1e-14)
				 U_V(m,n) = 0.0;
			printf("% 10.8f ", U_V(m,n));
		}
		cout << endl;
	}
	cout << endl;
}

int main()
{
	cout << "UNSYMMETRIC MATRIX:" << endl << endl;
	check(false);
	cout << endl;
	cout << "SYMMETRIC MATRIX:" << endl << endl;
	check(true);
}
