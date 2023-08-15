/*
 * moran_I.c
 *
 *  Created on: May 19, 2023
 *      Author: jiangp4
 */

#include "MoranI.h"

#include <stdio.h>
#include <stdlib.h>

// distance between spots
#define VISIUM_DISTANCE 100
#define OLD_DISTANCE 200
#define SQR(x) ((x) * (x))

double decay(const double d, const uint mode)
{
	if(mode == VISIUM){
		return exp(-SQR(d)/(2*SQR(75)));
	}else{
		return exp(-SQR(d)/(2*SQR(1.2)));
	}
}

void create_distance(double distance[], const uint max_shift, const uint mode)
{
	uint i, j;

	double x, y, d, visum_shift = 0.5*sqrt(3), old_shift = 0.5;

	for(i=0;i<max_shift; i++)
	{	// row shift

		for(j=0;j<2*max_shift;j++)
		{	// col shift

			if(mode == VISIUM){
				x = 0.5* j * VISIUM_DISTANCE;
				y = (i * visum_shift) * VISIUM_DISTANCE;

			}else if(mode == OLD){
				x = 0.5* j * OLD_DISTANCE;
				y = (i * old_shift) * OLD_DISTANCE;

			}else{
				fprintf(stderr, "Cannot recognize mode %u\n.", mode);
				exit(1);
			}

			d = sqrt(x*x + y*y);

			distance[distance_index(i, j, max_shift)] = decay(d,mode);
		}
	}
}

void z_normalize(double data[], const uint nrow, const uint ncol)
{
	uint i, j;
	double *ptr, aver, std;

	for (i=0, ptr = data; i< nrow; i++, ptr += ncol)
	{
		// compute average
		for(aver=0, j=0; j<ncol; j++){
			aver += ptr[j];
		}

		// compute standard deviation
		aver /= ncol;

		for(std=0, j=0; j<ncol; j++){
			std += (ptr[j] - aver)*(ptr[j] - aver);
		}

		std = sqrt(std/ncol);

		for(j=0; j<ncol; j++){
			ptr[j] = (ptr[j] - aver)/std;
		}
	}
}

double moran_I(
	// input parameter group 1: data
	double x[], double y[], const uint n,

	// input parameter group 2: pre-computed distance
	double distance[], const uint max_shift,

	// input parameter group 3: spot index
	uint spot_row[],
	uint spot_col[],
	const uint max_spot_row,
	const uint max_spot_col,
	const uint spot_index_map[]
	)
{
	uint i, j, row_offset, inx_row, col_offset, inx_col;

	double w, weight_sum = 0, total_sum = 0;

	for(i = 0; i< n; i++)
	{
		for (inx_row=spot_row[i], row_offset=0; inx_row < MIN(spot_row[i] + max_shift, max_spot_row); inx_row++, row_offset++)
		{
			if(row_offset == 0){
				inx_col= spot_col[i];
			}else{
				inx_col= (spot_col[i] + 1 > 2*max_shift)? (spot_col[i] + 1 - 2*max_shift): 0;
			}

			// winding for the first start spot
			for (;inx_col < MIN(spot_col[i] + 2*max_shift, max_spot_col); inx_col++)
			{
				j = spot_index_map[spot_index(inx_row, inx_col, max_spot_col)];
				if(j != UINT_MAX) break;
			}

			for (;inx_col < MIN(spot_col[i] + 2*max_shift, max_spot_col); inx_col+=2)
			{
				j = spot_index_map[spot_index(inx_row, inx_col, max_spot_col)];

				if(j != UINT_MAX)
				{
					col_offset = abs_offset(spot_col[i], inx_col);
					w = distance[distance_index(row_offset, col_offset, max_shift)];

					total_sum += w * x[i] * y[j] + w * x[j] * y[i];
					weight_sum += w * 2;
				}
			}
		}

		/*
		for(j=i; j<n; j++)
		{
			row_offset = abs_offset(spot_row[i], spot_row[j]);
			col_offset = abs_offset(spot_col[i], spot_col[j]);

			// will not reduce again
			if(	row_offset >= max_shift) break;

			if(	col_offset >= 2*max_shift ) continue;

			w = distance[distance_index(row_offset, col_offset, max_shift)];

			total_sum += w * x[i] * y[j];
			weight_sum += w;
		}*/
	}

	if(weight_sum == 0){
		fprintf(stderr, "Warning: no pairwise distance is included.\n");
		return 0;
	}

	return total_sum / weight_sum;
}
