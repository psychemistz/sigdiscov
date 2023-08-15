/*
 * Computer Moran I between x and y
 * Peng Jiang: pengj@alumni.princeton.edu
 *
 *  Created on: May 19, 2023
 */

#ifndef MORAN_I_
#define MORAN_I_

#define VISIUM 0
#define OLD 1

#define MAX(x, y) ((x)>(y)?(x):(y))
#define MIN(x, y) ((x)<(y)?(x):(y))

// map spot ID ixj to data matrix column index
#define spot_index(i, j, max_index) (i)*(max_index) + (j)

// convert row shift, col shift to index in the distance array
#define distance_index(i, j, max_shift) (i)*2*(max_shift) + (j)

// absolute offset between x and y
#define abs_offset(x, y) ((x)>(y)?((x)-(y)):((y)-(x)))

#include <math.h>
#include <limits.h>
#include <assert.h>

typedef unsigned int uint;


// distance decay funtion
double decay(const double d, const uint mode);

// compute distance function of ST grid
void create_distance(double distance[], const uint max_shift, const uint mode);

void z_normalize(double data[], const uint nrow, const uint ncol);

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
	);

#endif /* MORAN_I_ */
