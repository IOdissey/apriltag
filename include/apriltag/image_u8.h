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

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "math_util.h"

namespace apriltag
{
	struct image_u8_t
	{
		int width;
		int height;
		uint8_t* buf;
	};

	image_u8_t* image_u8_create(const int width, const int height)
	{
		image_u8_t* im = (image_u8_t*)malloc(sizeof(image_u8_t));
		im->width = width;
		im->height = height;
		im->buf = (uint8_t*)malloc(height * width * sizeof(uint8_t));
		return im;
	}

	image_u8_t* image_u8_copy(const image_u8_t* in)
	{
		image_u8_t* copy = image_u8_create(in->width, in->height);
		memcpy(copy->buf, &in->buf, in->height * in->width * sizeof(uint8_t));
		return copy;
	}

	void image_u8_destroy(image_u8_t* im)
	{
		if (!im)
			return;
		free(im->buf);
		free(im);
	}

	// 1.5, 2, 3, 4, ... supported
	image_u8_t* image_u8_decimate(const image_u8_t* im, float ffactor)
	{
		int width = im->width, height = im->height;

		if (ffactor == 1.5) {
			int swidth = width / 3 * 2, sheight = height / 3 * 2;

			image_u8_t* decim = image_u8_create(swidth, sheight);

			int y = 0, sy = 0;
			while (sy < sheight) {
				int x = 0, sx = 0;
				while (sx < swidth) {

					// a b c
					// d e f
					// g h i
					uint8_t a = im->buf[(y + 0) * im->width + (x + 0)];
					uint8_t b = im->buf[(y + 0) * im->width + (x + 1)];
					uint8_t c = im->buf[(y + 0) * im->width + (x + 2)];

					uint8_t d = im->buf[(y + 1) * im->width + (x + 0)];
					uint8_t e = im->buf[(y + 1) * im->width + (x + 1)];
					uint8_t f = im->buf[(y + 1) * im->width + (x + 2)];

					uint8_t g = im->buf[(y + 2) * im->width + (x + 0)];
					uint8_t h = im->buf[(y + 2) * im->width + (x + 1)];
					uint8_t i = im->buf[(y + 2) * im->width + (x + 2)];

					decim->buf[(sy + 0) * decim->width + (sx + 0)] =
						(4 * a + 2 * b + 2 * d + e) / 9;
					decim->buf[(sy + 0) * decim->width + (sx + 1)] =
						(4 * c + 2 * b + 2 * f + e) / 9;

					decim->buf[(sy + 1) * decim->width + (sx + 0)] =
						(4 * g + 2 * d + 2 * h + e) / 9;
					decim->buf[(sy + 1) * decim->width + (sx + 1)] =
						(4 * i + 2 * f + 2 * h + e) / 9;

					x += 3;
					sx += 2;
				}

				y += 3;
				sy += 2;
			}

			return decim;
		}

		int factor = (int)ffactor;

		int swidth = 1 + (width - 1) / factor;
		int sheight = 1 + (height - 1) / factor;
		image_u8_t* decim = image_u8_create(swidth, sheight);
		int sy = 0;
		for (int y = 0; y < height; y += factor) {
			int sx = 0;
			for (int x = 0; x < width; x += factor) {
				decim->buf[sy * decim->width + sx] = im->buf[y * im->width + x];
				sx++;
			}
			sy++;
		}
		return decim;
	}

	void convolve(const uint8_t* x, uint8_t* y, const int sz, const uint8_t* k, const int ksz)
	{
		assert((ksz & 1) == 1);

		const int ksz_2 = ksz / 2;
		for (int i = 0; i < ksz_2 && i < sz; i++)
			y[i] = x[i];

		for (int i = 0; i < sz - ksz; i++)
		{
			uint32_t acc = 0;
			for (int j = 0; j < ksz; j++)
				acc += k[j] * x[i + j];
			y[ksz_2 + i] = acc >> 8;
		}

		for (int i = sz - ksz + ksz_2; i < sz; i++)
			y[i] = x[i];
	}

	void image_u8_convolve_2D(image_u8_t* im, const uint8_t* k, const int ksz)
	{
		assert((ksz & 1) == 1); // ksz must be odd.

		const int h = im->height;
		const int w = im->width;
		const int wh = w * h;
		uint8_t* b1 = (uint8_t*)malloc(sizeof(uint8_t) * imax(w, 2 * h));
		uint8_t* b2 = b1 + h;

		for (int i = 0; i < wh; i += w)
		{
			memcpy(b1, &im->buf[i], w);
			convolve(b1, &im->buf[i], w, k, ksz);
		}

		for (int x = 0; x < w; x++)
		{
			for (int y = 0, i = x; y < h; y++, i += w)
				b1[y] = im->buf[i];
			convolve(b1, b2, h, k, ksz);
			for (int y = 0, i = x; y < h; y++, i += w)
				im->buf[i] = b2[y];
		}

		free(b1);
	}

	void image_u8_gaussian_blur(image_u8_t* im, double sigma, const int ksz)
	{
		if (sigma == 0)
			return;

		assert((ksz & 1) == 1); // ksz must be odd.

		// build the kernel.
		double* dk = (double*)malloc(sizeof(double) * ksz);

		// for kernel of length 5:
		// dk[0] = f(-2), dk[1] = f(-1), dk[2] = f(0), dk[3] = f(1), dk[4] = f(2)
		for (int i = 0; i < ksz; i++) {
			int x = -ksz / 2 + i;
			double v = exp(-.5 * sq(x / sigma));
			dk[i] = v;
		}

		// normalize
		double acc = 0;
		for (int i = 0; i < ksz; i++)
			acc += dk[i];

		for (int i = 0; i < ksz; i++)
			dk[i] /= acc;

		uint8_t* k = (uint8_t*)malloc(sizeof(uint8_t) * ksz);
		for (int i = 0; i < ksz; i++)
			k[i] = dk[i] * 255;

		free(dk);

		image_u8_convolve_2D(im, k, ksz);
		free(k);
	}
}