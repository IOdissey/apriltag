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

#include "zarray.h"
#include "math.h"


namespace apriltag
{
	////////////////////////////////////////////////////////////////////
	// Lines
	struct g2d_line_t
	{
		// Internal representation: a point that the line goes through (p) and
		// the direction of the line (u).
		double p[2];
		double u[2]; // always a unit vector
	};

	////////////////////////////////////////////////////////////////////
	// Line Segments. line.p is always one endpoint; p1 is the other
	// endpoint.
	struct g2d_line_segment_t
	{
		g2d_line_t line;
		double p1[2];
	};

	zarray_t* g2d_polygon_create_zeros(int sz)
	{
		zarray_t* points = zarray_create(sizeof(double[2]));

		double z[2] = { 0, 0 };

		for (int i = 0; i < sz; i++)
			zarray_add(points, z);

		return points;
	}

	void g2d_line_init_from_points(g2d_line_t* line, const double p0[2], const double p1[2])
	{
		line->p[0] = p0[0];
		line->p[1] = p0[1];
		line->u[0] = p1[0] - p0[0];
		line->u[1] = p1[1] - p0[1];
		double mag = sqrtf(sq(line->u[0]) + sq(line->u[1]));

		line->u[0] /= mag;
		line->u[1] /= mag;
	}

	void g2d_line_segment_init_from_points(g2d_line_segment_t* seg, const double p0[2], const double p1[2])
	{
		g2d_line_init_from_points(&seg->line, p0, p1);
		seg->p1[0] = p1[0];
		seg->p1[1] = p1[1];
	}

	// The line defines a one-dimensional coordinate system whose origin
	// is p. Where is q? (If q is not on the line, the point nearest q is
	// returned.
	double g2d_line_get_coordinate(const g2d_line_t* line, const double q[2])
	{
		return (q[0] - line->p[0]) * line->u[0] + (q[1] - line->p[1]) * line->u[1];
	}

	// Compute intersection of two line segments. If they intersect,
	// result is stored in p and 1 is returned. Otherwise, zero is
	// returned. p may be NULL.
	int g2d_line_intersect_line(const g2d_line_t* linea, const g2d_line_t* lineb, double* p)
	{
		// this implementation is many times faster than the original,
		// mostly due to avoiding a general-purpose LU decomposition in
		// Matrix.inverse().
		double m00, m01, m10, m11;
		double i00, i01;
		double b00, b10;

		m00 = linea->u[0];
		m01 = -lineb->u[0];
		m10 = linea->u[1];
		m11 = -lineb->u[1];

		// determinant of m
		double det = m00 * m11 - m01 * m10;

		// parallel lines?
		if (fabs(det) < 0.00000001)
			return 0;

		// inverse of m
		i00 = m11 / det;
		i01 = -m01 / det;

		b00 = lineb->p[0] - linea->p[0];
		b10 = lineb->p[1] - linea->p[1];

		double x00; //, x10;
		x00 = i00 * b00 + i01 * b10;

		if (p != NULL) {
			p[0] = linea->u[0] * x00 + linea->p[0];
			p[1] = linea->u[1] * x00 + linea->p[1];
		}

		return 1;
	}

	// Compute intersection of two line segments. If they intersect,
	// result is stored in p and 1 is returned. Otherwise, zero is
	// returned. p may be NULL.
	int g2d_line_segment_intersect_segment(const g2d_line_segment_t* sega, const g2d_line_segment_t* segb, double* p)
	{
		double tmp[2];

		if (!g2d_line_intersect_line(&sega->line, &segb->line, tmp))
			return 0;

		double a = g2d_line_get_coordinate(&sega->line, sega->line.p);
		double b = g2d_line_get_coordinate(&sega->line, sega->p1);
		double c = g2d_line_get_coordinate(&sega->line, tmp);

		// does intersection lie on the first line?
		if ((c < a && c < b) || (c > a && c > b))
			return 0;

		a = g2d_line_get_coordinate(&segb->line, segb->line.p);
		b = g2d_line_get_coordinate(&segb->line, segb->p1);
		c = g2d_line_get_coordinate(&segb->line, tmp);

		// does intersection lie on second line?
		if ((c < a && c < b) || (c > a && c > b))
			return 0;

		if (p != NULL) {
			p[0] = tmp[0];
			p[1] = tmp[1];
		}

		return 1;
	}

	// do the edges of polya and polyb collide? (Does NOT test for containment).
	int g2d_polygon_intersects_polygon(const zarray_t* polya, const zarray_t* polyb)
	{
		// do any of the line segments collide? If so, the answer is no.

		// dumb N^2 method.
		for (int ia = 0; ia < zarray_size(polya); ia++) {
			double pa0[2], pa1[2];
			zarray_get(polya, ia, pa0);
			zarray_get(polya, (ia + 1) % zarray_size(polya), pa1);

			g2d_line_segment_t sega;
			g2d_line_segment_init_from_points(&sega, pa0, pa1);

			for (int ib = 0; ib < zarray_size(polyb); ib++) {
				double pb0[2], pb1[2];
				zarray_get(polyb, ib, pb0);
				zarray_get(polyb, (ib + 1) % zarray_size(polyb), pb1);

				g2d_line_segment_t segb;
				g2d_line_segment_init_from_points(&segb, pb0, pb1);

				if (g2d_line_segment_intersect_segment(&sega, &segb, NULL))
					return 1;
			}
		}

		return 0;
	}

	// compute a point that is inside the polygon. (It may not be *far* inside though)
	void g2d_polygon_get_interior_point(const zarray_t* poly, double* p)
	{
		// take the first three points, which form a triangle. Find the middle point
		double a[2], b[2], c[2];

		zarray_get(poly, 0, a);
		zarray_get(poly, 1, b);
		zarray_get(poly, 2, c);

		p[0] = (a[0] + b[0] + c[0]) / 3;
		p[1] = (a[1] + b[1] + c[1]) / 3;
	}

	// Return 1 if point q lies within poly.
	int g2d_polygon_contains_point(const zarray_t* poly, double q[2])
	{
		// use winding. If the point is inside the polygon, we'll wrap
		// around it (accumulating 6.28 radians). If we're outside the
		// polygon, we'll accumulate zero.
		int psz = zarray_size(poly);
		assert(psz > 0);

		int last_quadrant;
		int quad_acc = 0;

		for (int i = 0; i <= psz; i++) {
			double* p;

			zarray_get_volatile(poly, i % psz, &p);

			// p[0] < q[0]       p[1] < q[1]    quadrant
			//     0                 0              0
			//     0                 1              3
			//     1                 0              1
			//     1                 1              2

			// p[1] < q[1]       p[0] < q[0]    quadrant
			//     0                 0              0
			//     0                 1              1
			//     1                 0              3
			//     1                 1              2

			int quadrant;
			if (p[0] < q[0])
				quadrant = (p[1] < q[1]) ? 2 : 1;
			else
				quadrant = (p[1] < q[1]) ? 3 : 0;

			if (i > 0) {
				int dquadrant = quadrant - last_quadrant;

				// encourage a jump table by mapping to small positive integers.
				switch (dquadrant) {
				case -3:
				case 1:
					quad_acc++;
					break;
				case -1:
				case 3:
					quad_acc--;
					break;
				case 0:
					break;
				case -2:
				case 2:
				{
					// get the previous point.
					double* p0;
					zarray_get_volatile(poly, i - 1, &p0);

					// Consider the points p0 and p (the points around the
					//polygon that we are tracing) and the query point q.
					//
					// If we've moved diagonally across quadrants, we want
					// to measure whether we have rotated +PI radians or
					// -PI radians. We can test this by computing the dot
					// product of vector (p0-q) with the vector
					// perpendicular to vector (p-q)
					double nx = p[1] - q[1];
					double ny = -p[0] + q[0];

					double dot = nx * (p0[0] - q[0]) + ny * (p0[1] - q[1]);
					if (dot < 0)
						quad_acc -= 2;
					else
						quad_acc += 2;

					break;
				}
				}
			}

			last_quadrant = quadrant;
		}

		int v = (quad_acc >= 2) || (quad_acc <= -2);

		return v;
	}

	// Is there some point which is in both polya and polyb?
	int g2d_polygon_overlaps_polygon(const zarray_t* polya, const zarray_t* polyb)
	{
		// do any of the line segments collide? If so, the answer is yes.
		if (g2d_polygon_intersects_polygon(polya, polyb))
			return 1;

		// if none of the edges cross, then the polygon is either fully
		// contained or fully outside.
		double p[2];
		g2d_polygon_get_interior_point(polyb, p);

		if (g2d_polygon_contains_point(polya, p))
			return 1;

		g2d_polygon_get_interior_point(polya, p);

		if (g2d_polygon_contains_point(polyb, p))
			return 1;

		return 0;
	}
}