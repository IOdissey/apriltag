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

#include <cstdint>

#include "zarray.h"
#include "image_u8.h"
#include "g2d.h"
#include "apriltag_types.h"
#include "apriltag_quad_thresh.h"
#include "math.h"


namespace apriltag
{
	void graymodel_init(graymodel_t* gm)
	{
		memset(gm, 0, sizeof(graymodel_t));
	}

	#define APRILTAG_U64_ONE ((uint64_t) 1)

	// don't forget to add a family!
	apriltag_detector_t* apriltag_detector_create()
	{
		apriltag_detector_t* td = (apriltag_detector_t*)calloc(1, sizeof(apriltag_detector_t));
		td->quad_decimate = 2.0;
		td->quad_sigma = 0.0;
		td->max_nmaxima = 10;
		td->min_cluster_pixels = 36;      // Some empirical value (> 0).
		td->max_cluster_pixels = 1000000; // Not limited.
		td->max_line_fit_mse = 10.0;
		td->cos_critical_rad = cos(20 * M_PI / 180);
		td->min_white_black_diff = 5;
		td->deglitch = 0;
		td->min_tag_area = 36;
		td->tag_families = zarray_create(sizeof(apriltag_family_t*));
		td->refine_edges = 1;
		td->decode_sharpening = 0.25;
		return td;
	}

	void quick_decode_add(struct quick_decode_t* qd, uint64_t code, int id, int hamming)
	{
		// uint32_t bucket = code % qd->nentries;
		uint32_t bucket = (code * 131) % qd->nentries;

		while (qd->entries[bucket].rcode != UINT64_MAX) {
			bucket = (bucket + 1) % qd->nentries;
		}

		qd->entries[bucket].rcode = code;
		qd->entries[bucket].id = id;
		qd->entries[bucket].hamming = hamming;
	}

	void quick_decode_init(apriltag_family_t* family, int maxhamming)
	{
		assert(family->impl == NULL);
		assert(family->ncodes < 65536);

		quick_decode_t* qd = (quick_decode_t*)calloc(1, sizeof(quick_decode_t));
		int capacity = family->ncodes;

		int nbits = family->nbits;

		if (maxhamming >= 1)
			capacity += family->ncodes * nbits;

		if (maxhamming >= 2)
			capacity += family->ncodes * nbits * (nbits - 1) / 2;

		if (maxhamming >= 3)
			capacity += family->ncodes * nbits * (nbits - 1) * (nbits - 2) / 6;

		qd->nentries = capacity * 3;

		qd->entries = (quick_decode_entry_t*)calloc(qd->nentries, sizeof(quick_decode_entry_t));
		if (qd->entries == NULL) {
			printf("apriltag.c: failed to allocate hamming decode table. Reduce max hamming size.\n");
			exit(-1);
		}

		for (int i = 0; i < qd->nentries; i++)
			qd->entries[i].rcode = UINT64_MAX;

		for (uint32_t i = 0; i < family->ncodes; i++) {
			uint64_t code = family->codes[i];

			// add exact code (hamming = 0)
			quick_decode_add(qd, code, i, 0);

			if (maxhamming >= 1) {
				// add hamming 1
				for (int j = 0; j < nbits; j++)
					quick_decode_add(qd, code ^ (APRILTAG_U64_ONE << j), i, 1);
			}

			if (maxhamming >= 2) {
				// add hamming 2
				for (int j = 0; j < nbits; j++)
					for (int k = 0; k < j; k++)
						quick_decode_add(qd, code ^ (APRILTAG_U64_ONE << j) ^ (APRILTAG_U64_ONE << k), i, 2);
			}

			if (maxhamming >= 3) {
				// add hamming 3
				for (int j = 0; j < nbits; j++)
					for (int k = 0; k < j; k++)
						for (int m = 0; m < k; m++)
							quick_decode_add(qd, code ^ (APRILTAG_U64_ONE << j) ^ (APRILTAG_U64_ONE << k) ^ (APRILTAG_U64_ONE << m), i, 3);
			}

			if (maxhamming > 3) {
				printf("apriltag.c: maxhamming beyond 3 not supported\n");
			}
		}

		family->impl = qd;
	}


	// add a family to the apriltag detector. caller still "owns" the family.
	// a single instance should only be provided to one apriltag detector instance.
	void apriltag_detector_add_family_bits(apriltag_detector_t* td, apriltag_family_t* fam, int bits_corrected)
	{
		zarray_add(td->tag_families, &fam);

		if (!fam->impl)
			quick_decode_init(fam, bits_corrected);
	}

	// Tunable, but really, 2 is a good choice. Values of >=3
	// consume prohibitively large amounts of memory, and otherwise
	// you want the largest value possible.
	static inline void apriltag_detector_add_family(apriltag_detector_t* td, apriltag_family_t* fam)
	{
		apriltag_detector_add_family_bits(td, fam, 2);
	}

	void quick_decode_uninit(apriltag_family_t* fam)
	{
		if (!fam->impl)
			return;

		quick_decode_t* qd = (quick_decode_t*) fam->impl;
		free(qd->entries);
		free(qd);
		fam->impl = NULL;
	}

	// does not deallocate the family.
	void apriltag_detector_remove_family(apriltag_detector_t* td, apriltag_family_t* fam)
	{
		quick_decode_uninit(fam);
		zarray_remove_value(td->tag_families, &fam, 0);
	}

	// unregister all families, but does not deallocate the underlying tag family objects.
	void apriltag_detector_clear_families(apriltag_detector_t* td)
	{
		for (int i = 0; i < zarray_size(td->tag_families); i++) {
			apriltag_family_t* fam;
			zarray_get(td->tag_families, i, &fam);
			quick_decode_uninit(fam);
		}
		zarray_clear(td->tag_families);
	}

	// Destroy the april tag detector (but not the underlying
	// apriltag_family_t used to initialize it.)
	void apriltag_detector_destroy(apriltag_detector_t* td)
	{
		apriltag_detector_clear_families(td);
		zarray_destroy(td->tag_families);
		free(td);
	}

	void refine_edges(const apriltag_detector_t* td, const image_u8_t* im_orig, quad_t* quad)
	{
		double lines[4][4]; // for each line, [Ex Ey nx ny]

		for (int edge = 0; edge < 4; edge++) {
			int a = edge, b = (edge + 1) & 3; // indices of the end points.

			// compute the normal to the current line estimate
			double nx = quad->p[b][1] - quad->p[a][1];
			double ny = -quad->p[b][0] + quad->p[a][0];
			double mag = sqrt(nx * nx + ny * ny);
			nx /= mag;
			ny /= mag;

			if (quad->reversed_border) {
				nx = -nx;
				ny = -ny;
			}

			// we will now fit a NEW line by sampling points near
			// our original line that have large gradients. On really big tags,
			// we're willing to sample more to get an even better estimate.
			int nsamples = imax(16, mag / 8); // XXX tunable

			// stats for fitting a line...
			double Mx = 0, My = 0, Mxx = 0, Mxy = 0, Myy = 0, N = 0;

			for (int s = 0; s < nsamples; s++) {
				// compute a point along the line... Note, we're avoiding
				// sampling *right* at the corners, since those points are
				// the least reliable.
				double alpha = (1.0 + s) / (nsamples + 1);
				double x0 = alpha * quad->p[a][0] + (1 - alpha) * quad->p[b][0];
				double y0 = alpha * quad->p[a][1] + (1 - alpha) * quad->p[b][1];

				// search along the normal to this line, looking at the
				// gradients along the way. We're looking for a strong
				// response.
				double Mn = 0;
				double Mcount = 0;

				// XXX tunable: how far to search?  We want to search far
				// enough that we find the best edge, but not so far that
				// we hit other edges that aren't part of the tag. We
				// shouldn't ever have to search more than quad_decimate,
				// since otherwise we would (ideally) have started our
				// search on another pixel in the first place. Likewise,
				// for very small tags, we don't want the range to be too
				// big.
				double range = td->quad_decimate + 1;

				// XXX tunable step size.
				for (double n = -range; n <= range; n += 0.25) {
					// Because of the guaranteed winding order of the
					// points in the quad, we will start inside the white
					// portion of the quad and work our way outward.
					//
					// sample to points (x1,y1) and (x2,y2) XXX tunable:
					// how far +/- to look? Small values compute the
					// gradient more precisely, but are more sensitive to
					// noise.
					double grange = 1;
					int x1 = x0 + (n + grange) * nx;
					int y1 = y0 + (n + grange) * ny;
					if (x1 < 0 || x1 >= im_orig->width || y1 < 0 || y1 >= im_orig->height)
						continue;

					int x2 = x0 + (n - grange) * nx;
					int y2 = y0 + (n - grange) * ny;
					if (x2 < 0 || x2 >= im_orig->width || y2 < 0 || y2 >= im_orig->height)
						continue;

					int g1 = im_orig->buf[y1 * im_orig->width + x1];
					int g2 = im_orig->buf[y2 * im_orig->width + x2];

					if (g1 < g2) // reject points whose gradient is "backwards". They can only hurt us.
						continue;

					double weight = (g2 - g1) * (g2 - g1); // XXX tunable. What shape for weight=f(g2-g1)?

					// compute weighted average of the gradient at this point.
					Mn += weight * n;
					Mcount += weight;
				}

				// what was the average point along the line?
				if (Mcount == 0)
					continue;

				double n0 = Mn / Mcount;

				// where is the point along the line?
				double bestx = x0 + n0 * nx;
				double besty = y0 + n0 * ny;

				// update our line fit statistics
				Mx += bestx;
				My += besty;
				Mxx += bestx * bestx;
				Mxy += bestx * besty;
				Myy += besty * besty;
				N++;
			}

			// fit a line
			double Ex = Mx / N, Ey = My / N;
			double Cxx = Mxx / N - Ex * Ex;
			double Cxy = Mxy / N - Ex * Ey;
			double Cyy = Myy / N - Ey * Ey;

			// TODO: Can replace this with same code as in fit_line.
			double normal_theta = .5 * atan2f(-2 * Cxy, (Cyy - Cxx));
			nx = cosf(normal_theta);
			ny = sinf(normal_theta);
			lines[edge][0] = Ex;
			lines[edge][1] = Ey;
			lines[edge][2] = nx;
			lines[edge][3] = ny;
		}

		// now refit the corners of the quad
		for (int i = 0; i < 4; i++) {

			const int i2 = (i + 1) & 3;

			// solve for the intersection of lines (i) and (i+1)&3.
			double A00 = lines[i][3], A01 = -lines[i2][3];
			double A10 = -lines[i][2], A11 = lines[i2][2];
			double B0 = -lines[i][0] + lines[i2][0];
			double B1 = -lines[i][1] + lines[i2][1];

			double det = A00 * A11 - A10 * A01;

			// inverse.
			if (fabs(det) > 0.001) {
				// solve
				double W00 = A11 / det, W01 = -A01 / det;

				double L0 = W00 * B0 + W01 * B1;

				// Compute intersection. Note that line i represents the line from corner i to (i+1)&3, so
				// the intersection of line i with line (i+1)&3 represents corner (i+1)&3.
				quad->p[i2][0] = lines[i][0] + L0 * A00;
				quad->p[i2][1] = lines[i][1] + L0 * A10;
			}
		}
	}

	// returns non-zero if an error occurs (i.e., H has no inverse)
	int quad_update_homographies(quad_t* quad)
	{
		double corr_arr[4][4];

		for (int i = 0; i < 4; i++) {
			corr_arr[i][0] = (i == 0 || i == 3) ? -1 : 1;
			corr_arr[i][1] = (i == 0 || i == 1) ? -1 : 1;
			corr_arr[i][2] = quad->p[i][0];
			corr_arr[i][3] = quad->p[i][1];
		}

		if (homography_compute(corr_arr, quad->h))
			return 0; // Success!

		return -1;
	}

	void graymodel_add(graymodel_t* gm, double x, double y, double gray)
	{
		// update upper right entries of A = J'J
		gm->A[0][0] += x * x;
		gm->A[0][1] += x * y;
		gm->A[0][2] += x;
		gm->A[1][1] += y * y;
		gm->A[1][2] += y;
		gm->A[2][2] += 1;

		// update B = J'gray
		gm->B[0] += x * gray;
		gm->B[1] += y * gray;
		gm->B[2] += gray;
	}

	void graymodel_solve(graymodel_t* gm)
	{
		mat33_sym_solve((double*)gm->A, gm->B, gm->C);
	}

	double graymodel_interpolate(const graymodel_t* gm, double x, double y)
	{
		return gm->C[0] * x + gm->C[1] * y + gm->C[2];
	}

	double value_for_pixel(const image_u8_t* im, double px, double py)
	{
		int x1 = floor(px - 0.5);
		int x2 = ceil(px - 0.5);
		double x = px - 0.5 - x1;
		int y1 = floor(py - 0.5);
		int y2 = ceil(py - 0.5);
		double y = py - 0.5 - y1;
		if (x1 < 0 || x2 >= im->width || y1 < 0 || y2 >= im->height) {
			return -1;
		}
		return im->buf[y1 * im->width + x1] * (1 - x) * (1 - y) +
			im->buf[y1 * im->width + x2] * x * (1 - y) +
			im->buf[y2 * im->width + x1] * (1 - x) * y +
			im->buf[y2 * im->width + x2] * x * y;
	}

	void sharpen(const apriltag_detector_t* td, double* values, int size)
	{
		double* sharpened = (double*)malloc(sizeof(double) * size * size);
		double kernel[9] = {
			0, -1, 0,
			-1, 4, -1,
			0, -1, 0
		};

		for (int y = 0; y < size; y++) {
			for (int x = 0; x < size; x++) {
				sharpened[y * size + x] = 0;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						if ((y + i - 1) < 0 || (y + i - 1) > size - 1 || (x + j - 1) < 0 || (x + j - 1) > size - 1) {
							continue;
						}
						sharpened[y * size + x] += values[(y + i - 1) * size + (x + j - 1)] * kernel[i * 3 + j];
					}
				}
			}
		}


		for (int y = 0; y < size; y++) {
			for (int x = 0; x < size; x++) {
				values[y * size + x] = values[y * size + x] + td->decode_sharpening * sharpened[y * size + x];
			}
		}

		free(sharpened);
	}

	/**
	 * Assuming we are drawing the image one quadrant at a time, what would the rotated image look like?
	 * Special care is taken to handle the case where there is a middle pixel of the image.
	 */
	uint64_t rotate90(uint64_t w, int numBits)
	{
		int p = numBits;
		uint64_t l = 0;
		if (numBits % 4 == 1) {
			p = numBits - 1;
			l = 1;
		}
		w = ((w >> l) << (p / 4 + l)) | (w >> (3 * p / 4 + l) << l) | (w & l);
		w &= ((APRILTAG_U64_ONE << numBits) - 1);
		return w;
	}

	// returns an entry with hamming set to 255 if no decode was found.
	void quick_decode_codeword(const apriltag_family_t* tf, uint64_t rcode, quick_decode_entry_t* entry)
	{
		quick_decode_t* qd = (quick_decode_t*)tf->impl;

		for (int ridx = 0; ridx < 4; ridx++) {

			for (int bucket = rcode % qd->nentries;
				qd->entries[bucket].rcode != UINT64_MAX;
				bucket = (bucket + 1) % qd->nentries) {

				if (qd->entries[bucket].rcode == rcode) {
					*entry = qd->entries[bucket];
					entry->rotation = ridx;
					return;
				}
			}

			rcode = rotate90(rcode, tf->nbits);
		}

		entry->rcode = 0;
		entry->id = 65535;
		entry->hamming = 255;
		entry->rotation = 0;
	}

	// returns the decision margin. Return < 0 if the detection should be rejected.
	float quad_decode(const apriltag_detector_t* td, const apriltag_family_t* family, const image_u8_t* im, const quad_t* quad, quick_decode_entry_t* entry)
	{
		// decode the tag binary contents by sampling the pixel
		// closest to the center of each bit cell.

		// We will compute a threshold by sampling known white/black cells around this tag.
		// This sampling is achieved by considering a set of samples along lines.
		//
		// coordinates are given in bit coordinates. ([0, fam->border_width]).
		//
		// { initial x, initial y, delta x, delta y, WHITE=1 }
		float patterns[] = {
			// left white column
			-0.5, 0.5,
			0, 1,
			1,

			// left black column
			0.5, 0.5,
			0, 1,
			0,

			// right white column
			family->width_at_border + 0.5f, .5,
			0, 1,
			1,

			// right black column
			family->width_at_border - 0.5f, .5,
			0, 1,
			0,

			// top white row
			0.5, -0.5,
			1, 0,
			1,

			// top black row
			0.5, 0.5,
			1, 0,
			0,

			// bottom white row
			0.5, family->width_at_border + 0.5f,
			1, 0,
			1,

			// bottom black row
			0.5, family->width_at_border - 0.5f,
			1, 0,
			0

			// XXX double-counts the corners.
		};

		graymodel_t whitemodel, blackmodel;
		graymodel_init(&whitemodel);
		graymodel_init(&blackmodel);

		for (int pattern_idx = 0; pattern_idx < sizeof(patterns) / (5 * sizeof(float)); pattern_idx++) {
			float* pattern = &patterns[pattern_idx * 5];

			int is_white = pattern[4];

			for (int i = 0; i < family->width_at_border; i++) {
				double tagx01 = (pattern[0] + i * pattern[2]) / (family->width_at_border);
				double tagy01 = (pattern[1] + i * pattern[3]) / (family->width_at_border);

				double tagx = 2 * (tagx01 - 0.5);
				double tagy = 2 * (tagy01 - 0.5);

				double px, py;
				homography_project(quad->h, tagx, tagy, &px, &py);

				// don't round
				int ix = px;
				int iy = py;
				if (ix < 0 || iy < 0 || ix >= im->width || iy >= im->height)
					continue;

				int v = im->buf[iy * im->width + ix];

				if (is_white)
					graymodel_add(&whitemodel, tagx, tagy, v);
				else
					graymodel_add(&blackmodel, tagx, tagy, v);
			}
		}

		if (family->width_at_border > 1) {
			graymodel_solve(&whitemodel);
			graymodel_solve(&blackmodel);
		}
		else {
			graymodel_solve(&whitemodel);
			blackmodel.C[0] = 0;
			blackmodel.C[1] = 0;
			blackmodel.C[2] = blackmodel.B[2] / 4;
		}

		// XXX Tunable
		if ((graymodel_interpolate(&whitemodel, 0, 0) - graymodel_interpolate(&blackmodel, 0, 0) < 0) != family->reversed_border) {
			return -1;
		}

		// compute the average decision margin (how far was each bit from
		// the decision boundary?
		//
		// we score this separately for white and black pixels and return
		// the minimum average threshold for black/white pixels. This is
		// to penalize thresholds that are too close to an extreme.
		float black_score = 0, white_score = 0;
		float black_score_count = 1, white_score_count = 1;

		double* values = (double*)calloc(family->total_width * family->total_width, sizeof(double));

		int min_coord = (family->width_at_border - family->total_width) / 2;
		for (uint32_t i = 0; i < family->nbits; i++) {
			int bity = family->bit_y[i];
			int bitx = family->bit_x[i];

			double tagx01 = (bitx + 0.5) / (family->width_at_border);
			double tagy01 = (bity + 0.5) / (family->width_at_border);

			// scale to [-1, 1]
			double tagx = 2 * (tagx01 - 0.5);
			double tagy = 2 * (tagy01 - 0.5);

			double px, py;
			homography_project(quad->h, tagx, tagy, &px, &py);

			double v = value_for_pixel(im, px, py);

			if (v == -1) {
				continue;
			}

			double thresh = (graymodel_interpolate(&blackmodel, tagx, tagy) + graymodel_interpolate(&whitemodel, tagx, tagy)) / 2.0;
			values[family->total_width * (bity - min_coord) + bitx - min_coord] = v - thresh;
		}

		sharpen(td, values, family->total_width);

		uint64_t rcode = 0;
		for (uint32_t i = 0; i < family->nbits; i++) {
			int bity = family->bit_y[i];
			int bitx = family->bit_x[i];
			rcode = (rcode << 1);
			double v = values[(bity - min_coord) * family->total_width + bitx - min_coord];

			if (v > 0) {
				white_score += v;
				white_score_count++;
				rcode |= 1;
			}
			else {
				black_score -= v;
				black_score_count++;
			}
		}

		quick_decode_codeword(family, rcode, entry);
		free(values);
		return fmin(white_score / white_score_count, black_score / black_score_count);
	}

	int prefer_smaller(int pref, double q0, double q1)
	{
		if (pref)     // already prefer something? exit.
			return pref;

		if (q0 < q1)
			return -1; // we now prefer q0
		if (q1 < q0)
			return 1; // we now prefer q1

		// no preference
		return 0;
	}

	inline int detection_compare_function(const void* _a, const void* _b)
	{
		apriltag_detection_t* a = *(apriltag_detection_t**)_a;
		apriltag_detection_t* b = *(apriltag_detection_t**)_b;

		return a->id - b->id;
	}

	// Detect tags from an image and return an array of
	// apriltag_detection_t*. You can use apriltag_detections_destroy to
	// free the array and the detections it contains, or call
	// _detection_destroy and zarray_destroy yourself.
	zarray_t* apriltag_detector_detect(apriltag_detector_t* td, image_u8_t* im_orig)
	{
		if (zarray_size(td->tag_families) == 0) {
			zarray_t* s = zarray_create(sizeof(apriltag_detection_t*));
			printf("No tag families enabled.");
			return s;
		}

		///////////////////////////////////////////////////////////
		// Step 1. Detect quads according to requested image decimation
		// and blurring parameters.
		image_u8_t* quad_im = im_orig;
		if (td->quad_decimate > 1) {
			quad_im = image_u8_decimate(im_orig, td->quad_decimate);
		}

		if (td->quad_sigma != 0) {
			// compute a reasonable kernel width by figuring that the
			// kernel should go out 2 std devs.
			//
			// max sigma          ksz
			// 0.499              1  (disabled)
			// 0.999              3
			// 1.499              5
			// 1.999              7

			float sigma = fabsf((float)td->quad_sigma);

			int ksz = 4 * sigma; // 2 std devs in each direction
			if ((ksz & 1) == 0)
				ksz++;

			if (ksz > 1) {

				if (td->quad_sigma > 0) {
					// Apply a blur
					image_u8_gaussian_blur(quad_im, sigma, ksz);
				}
				else {
					// SHARPEN the image by subtracting the low frequency components.
					image_u8_t* orig = image_u8_copy(quad_im);
					image_u8_gaussian_blur(quad_im, sigma, ksz);

					for (int y = 0; y < orig->height; y++) {
						for (int x = 0; x < orig->width; x++) {
							int vorig = orig->buf[y * orig->width + x];
							int vblur = quad_im->buf[y * quad_im->width + x];

							int v = 2 * vorig - vblur;
							if (v < 0)
								v = 0;
							if (v > 255)
								v = 255;

							quad_im->buf[y * quad_im->width + x] = (uint8_t)v;
						}
					}
					image_u8_destroy(orig);
				}
			}
		}

		zarray_t* quads = apriltag_quad_thresh(td, quad_im);

		// adjust centers of pixels so that they correspond to the
		// original full-resolution image.
		if (td->quad_decimate > 1) {
			for (int i = 0; i < zarray_size(quads); i++) {
				quad_t* q;
				zarray_get_volatile(quads, i, &q);

				for (int j = 0; j < 4; j++) {
					if (td->quad_decimate == 1.5) {
						q->p[j][0] *= td->quad_decimate;
						q->p[j][1] *= td->quad_decimate;
					}
					else {
						q->p[j][0] = (q->p[j][0] - 0.5) * td->quad_decimate + 0.5;
						q->p[j][1] = (q->p[j][1] - 0.5) * td->quad_decimate + 0.5;
					}
				}
			}
		}

		if (quad_im != im_orig)
			image_u8_destroy(quad_im);

		zarray_t* detections = zarray_create(sizeof(apriltag_detection_t*));

		td->nquads = zarray_size(quads);

		////////////////////////////////////////////////////////////////
		// Step 2. Decode tags from each quad.
		const int quadidx_max = zarray_size(quads);
		for (int quadidx = 0; quadidx < quadidx_max; quadidx++)
		{
			quad_t* quad_original;
			zarray_get_volatile(quads, quadidx, &quad_original);

			// refine edges is not dependent upon the tag family, thus
			// apply this optimization BEFORE the other work.
			//if (td->quad_decimate > 1 && td->refine_edges) {
			if (td->refine_edges)
				refine_edges(td, im_orig, quad_original);

			// make sure the homographies are computed...
			if (quad_update_homographies(quad_original) != 0)
				continue;

			for (int famidx = 0; famidx < zarray_size(td->tag_families); famidx++)
			{
				apriltag_family_t* family;
				zarray_get(td->tag_families, famidx, &family);

				if (family->reversed_border != quad_original->reversed_border)
					continue;

				// since the geometry of tag families can vary, start any
				// optimization process over with the original quad.
				quad_t* quad = quad_copy(quad_original);

				quick_decode_entry_t entry;
				float decision_margin = quad_decode(td, family, im_orig, quad, &entry);

				if (decision_margin >= 0 && entry.hamming < 255)
				{
					apriltag_detection_t* det = (apriltag_detection_t*)calloc(1, sizeof(apriltag_detection_t));

					det->family = family;
					det->id = entry.id;
					det->hamming = entry.hamming;
					det->decision_margin = decision_margin;

					double theta = entry.rotation * M_PI / 2.0;
					double c = cos(theta), s = sin(theta);

					// Fix the rotation of our homography to properly orient the tag
					double r[] = {
						c, -s,  0,
						s,  c,  0,
						0,  0,  1
					};
					mat33_mull(quad->h, r, det->h);

					homography_project(det->h, 0, 0, &det->c[0], &det->c[1]);

					// [-1, -1], [1, -1], [1, 1], [-1, 1], Desired points
					// [-1, 1], [1, 1], [1, -1], [-1, -1], FLIP Y
					// adjust the points in det->p so that they correspond to
					// counter-clockwise around the quad, starting at -1,-1.
					for (int i = 0; i < 4; i++)
					{
						int tcx = (i == 1 || i == 2) ? 1 : -1;
						int tcy = (i < 2) ? 1 : -1;

						double p[2];

						// homography_project(det->H, tcx, tcy, &p[0], &p[1]);
						homography_project(det->h, tcx, tcy, &p[0], &p[1]);

						det->p[i][0] = p[0];
						det->p[i][1] = p[1];
					}
					zarray_add(detections, &det);
				}
				quad_destroy(quad);
			}
		}
		for (int i = 0; i < zarray_size(quads); i++)
		{
			quad_t* quad;
			zarray_get_volatile(quads, i, &quad);
		}
		zarray_destroy(quads);

		////////////////////////////////////////////////////////////////
		// Step 3. Reconcile detections--- don't report the same tag more
		// than once. (Allow non-overlapping duplicate detections.)
		zarray_t* poly0 = g2d_polygon_create_zeros(4);
		zarray_t* poly1 = g2d_polygon_create_zeros(4);
		for (int i0 = 0; i0 < zarray_size(detections); i0++)
		{
			apriltag_detection_t* det0;
			zarray_get(detections, i0, &det0);

			for (int k = 0; k < 4; k++)
				zarray_set(poly0, k, det0->p[k], NULL);

			for (int i1 = i0 + 1; i1 < zarray_size(detections); i1++)
			{
				apriltag_detection_t* det1;
				zarray_get(detections, i1, &det1);

				if (det0->id != det1->id || det0->family != det1->family)
					continue;

				for (int k = 0; k < 4; k++)
					zarray_set(poly1, k, det1->p[k], NULL);

				if (g2d_polygon_overlaps_polygon(poly0, poly1))
				{
					// the tags overlap. Delete one, keep the other.

					int pref = 0; // 0 means undecided which one we'll keep.
					pref = prefer_smaller(pref, det0->hamming, det1->hamming);     // want small hamming
					pref = prefer_smaller(pref, -det0->decision_margin, -det1->decision_margin);      // want bigger margins

					// if we STILL don't prefer one detection over the other, then pick
					// any deterministic criterion.
					for (int i = 0; i < 4; i++)
					{
						pref = prefer_smaller(pref, det0->p[i][0], det1->p[i][0]);
						pref = prefer_smaller(pref, det0->p[i][1], det1->p[i][1]);
					}

					if (pref == 0)
					{
						// at this point, we should only be undecided if the tag detections
						// are *exactly* the same. How would that happen?
						printf("uh oh, no preference for overlappingdetection\n");
					}

					if (pref < 0)
					{
						// keep det0, destroy det1
						apriltag_detection_destroy(det1);
						zarray_remove_index(detections, i1, 1);
						i1--; // retry the same index
						// goto retry1;
						continue;
					}
					else
					{
						// keep det1, destroy det0
						apriltag_detection_destroy(det0);
						zarray_remove_index(detections, i0, 1);
						i0--; // retry the same index.
						// goto retry0;
						break;
					}
				}
			// retry1:;
			}
		// retry0:;
		}
		zarray_destroy(poly0);
		zarray_destroy(poly1);

		zarray_sort(detections, detection_compare_function);

		return detections;
	}
}