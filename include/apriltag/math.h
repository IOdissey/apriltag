/*
Copyright (C) 2013-2016, The Regents of The University of Michigan.
All rights reserved.
This software was developed in the APRIL Robotics Lab under the
direction of Edwin Olson, ebolson@umich.edu. This software may be
available under alternative licensing terms; contact the address above.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the Regents of The University of Michigan.
*/

#pragma once

#include <cmath>


namespace apriltag
{
	#ifndef M_PI
	#define M_PI 3.141592653589793238462643383279502884196
	#endif
	#define M_TWOPI 6.2831853071795862319959

	inline double sq(double v)
	{
		return v * v;
	}

	inline int imin(int a, int b)
	{
		return (a < b) ? a : b;
	}

	inline int imax(int a, int b)
	{
		return (a > b) ? a : b;
	}

	// Inverse matrix (3x3).
	// b = a^-1
	//
	// | a0 a1 a2 |-1      1    | b0 b1 b2 |
	// | a3 a4 a5 |   = ------- | b3 b4 b5 |
	// | a6 a7 a8 |      det(A) | b6 b7 b8 |
	//
	bool mat33_inv(const double* a, double* b)
	{
		b[0] = a[4] * a[8] - a[5] * a[7];
		b[1] = a[2] * a[7] - a[1] * a[8];
		b[2] = a[1] * a[5] - a[2] * a[4];
		b[3] = a[5] * a[6] - a[3] * a[8];
		b[4] = a[0] * a[8] - a[2] * a[6];
		b[5] = a[2] * a[3] - a[0] * a[5];
		b[6] = a[3] * a[7] - a[4] * a[6];
		b[7] = a[1] * a[6] - a[0] * a[7];
		b[8] = a[0] * a[4] - a[1] * a[3];
		const double det = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
		// Determinant too small.
		if (std::abs(det) < 1e-9)
			return false;
		const double det_inv = 1.0 / det;
		for (int i = 0; i < 9; ++i)
			b[i] *= det_inv;
		return true;
	}

	// Matrix multiplication (3x3).
	// c = a * b
	//
	// | c0 c1 c2 |   | a0 a1 a2 |   | b0 b1 b2 |
	// | c3 c4 c5 | = | a3 a4 a5 | * | b3 b4 b5 |
	// | c6 c7 c8 |   | a6 a7 a8 |   | b6 b7 b8 |
	//
	void mat33_mull(const double* a, const double* b, double* c)
	{
		c[0] = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
		c[1] = a[0] * b[1] + a[1] * b[4] + a[2] * b[7];
		c[2] = a[0] * b[2] + a[1] * b[5] + a[2] * b[8];
		c[3] = a[3] * b[0] + a[4] * b[3] + a[5] * b[6];
		c[4] = a[3] * b[1] + a[4] * b[4] + a[5] * b[7];
		c[5] = a[3] * b[2] + a[4] * b[5] + a[5] * b[8];
		c[6] = a[6] * b[0] + a[7] * b[3] + a[8] * b[6];
		c[7] = a[6] * b[1] + a[7] * b[4] + a[8] * b[7];
		c[8] = a[6] * b[2] + a[7] * b[5] + a[8] * b[8];
	}

	// Computes the cholesky factorization of A, putting the lower
	// triangular matrix into R.
	inline void mat33_chol(const double* A, double* R)
	{
		// A[0] = R[0]*R[0]
		R[0] = sqrt(A[0]);

		// A[1] = R[0]*R[3];
		R[3] = A[1] / R[0];

		// A[2] = R[0]*R[6];
		R[6] = A[2] / R[0];

		// A[4] = R[3]*R[3] + R[4]*R[4]
		R[4] = sqrt(A[4] - R[3] * R[3]);

		// A[5] = R[3]*R[6] + R[4]*R[7]
		R[7] = (A[5] - R[3] * R[6]) / R[4];

		// A[8] = R[6]*R[6] + R[7]*R[7] + R[8]*R[8]
		R[8] = sqrt(A[8] - R[6] * R[6] - R[7] * R[7]);

		R[1] = 0;
		R[2] = 0;
		R[5] = 0;
	}

	inline void mat33_lower_tri_inv(const double* A, double* R)
	{
		// A[0]*R[0] = 1
		R[0] = 1 / A[0];

		// A[3]*R[0] + A[4]*R[3] = 0
		R[3] = -A[3] * R[0] / A[4];

		// A[4]*R[4] = 1
		R[4] = 1 / A[4];

		// A[6]*R[0] + A[7]*R[3] + A[8]*R[6] = 0
		R[6] = (-A[6] * R[0] - A[7] * R[3]) / A[8];

		// A[7]*R[4] + A[8]*R[7] = 0
		R[7] = -A[7] * R[4] / A[8];

		// A[8]*R[8] = 1
		R[8] = 1 / A[8];
	}

	inline void mat33_sym_solve(const double* A, const double* B, double* R)
	{
		double L[9];
		mat33_chol(A, L);

		double M[9];
		mat33_lower_tri_inv(L, M);

		double tmp[3];
		tmp[0] = M[0] * B[0];
		tmp[1] = M[3] * B[0] + M[4] * B[1];
		tmp[2] = M[6] * B[0] + M[7] * B[1] + M[8] * B[2];

		R[0] = M[0] * tmp[0] + M[3] * tmp[1] + M[6] * tmp[2];
		R[1] = M[4] * tmp[1] + M[7] * tmp[2];
		R[2] = M[8] * tmp[2];
	}

	//
	bool homography_compute(double c[4][4], double* h)
	{
		double a[] = {
			c[0][0], c[0][1], 1,       0,       0, 0, -c[0][0] * c[0][2], -c[0][1] * c[0][2], c[0][2],
			      0,       0, 0, c[0][0], c[0][1], 1, -c[0][0] * c[0][3], -c[0][1] * c[0][3], c[0][3],
			c[1][0], c[1][1], 1,       0,       0, 0, -c[1][0] * c[1][2], -c[1][1] * c[1][2], c[1][2],
			      0,       0, 0, c[1][0], c[1][1], 1, -c[1][0] * c[1][3], -c[1][1] * c[1][3], c[1][3],
			c[2][0], c[2][1], 1,       0,       0, 0, -c[2][0] * c[2][2], -c[2][1] * c[2][2], c[2][2],
			      0,       0, 0, c[2][0], c[2][1], 1, -c[2][0] * c[2][3], -c[2][1] * c[2][3], c[2][3],
			c[3][0], c[3][1], 1,       0,       0, 0, -c[3][0] * c[3][2], -c[3][1] * c[3][2], c[3][2],
			      0,       0, 0, c[3][0], c[3][1], 1, -c[3][0] * c[3][3], -c[3][1] * c[3][3], c[3][3],
		};
		// Eliminate.
		for (int col = 0; col < 8; col++)
		{
			// Find best row to swap with.
			double max_val = 0;
			int max_val_idx = -1;
			for (int row = col; row < 8; row++)
			{
				double val = std::abs(a[row * 9 + col]);
				if (val > max_val)
				{
					max_val = val;
					max_val_idx = row;
				}
			}
			if (max_val < 1e-9)
				return false;
			// Swap to get best row.
			if (max_val_idx != col)
			{
				for (int i = col; i < 9; i++)
				{
					double tmp = a[col * 9 + i];
					a[col * 9 + i] = a[max_val_idx * 9 + i];
					a[max_val_idx * 9 + i] = tmp;
				}
			}
			// Do eliminate.
			for (int i = col + 1; i < 8; i++)
			{
				double f = a[i * 9 + col] / a[col * 9 + col];
				a[i * 9 + col] = 0;
				for (int j = col + 1; j < 9; j++)
					a[i * 9 + j] -= f * a[col * 9 + j];
			}
		}
		// Back solve.
		for (int col = 7; col >= 0; col--)
		{
			double sum = 0;
			for (int i = col + 1; i < 8; i++)
				sum += a[col * 9 + i] * a[i * 9 + 8];
			a[col * 9 + 8] = (a[col * 9 + 8] - sum) / a[col * 9 + col];
		}
		h[0] = a[8];
		h[1] = a[17];
		h[2] = a[26];
		h[3] = a[35];
		h[4] = a[44];
		h[5] = a[53];
		h[6] = a[62];
		h[7] = a[71];
		h[8] = 1.0;
		return true;
	}

	// h - 3x3 dimansion.
	void homography_project(const double* h, double x, double y, double *ox, double *oy)
	{
		double xx = h[0] * x + h[1] * y + h[2];
		double yy = h[3] * x + h[4] * y + h[5];
		double zz = h[6] * x + h[7] * y + h[8];
		*ox = xx / zz;
		*oy = yy / zz;
	}
}