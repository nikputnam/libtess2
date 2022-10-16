#include "bez.h"
void cubicBez(float t, float* p0, float* p1, float* p2, float* p3, float* x) {
	// B(t) = (1-t)^3 P0 + 3*(1-t)^2 t P1 + 3*(1-t)*t^2 P2 + t^3 P3
	//   
	x[0]=0.0;
	x[1]=0.0;

	float t_squared = t*t;
	float t_cubed = t*t_squared;
	float one_minus_t = 1.0-t;
	float one_minus_t_squared = one_minus_t*one_minus_t;
	float one_minus_t_cubed = one_minus_t*one_minus_t_squared;

	x[0] += one_minus_t_cubed * p0[0];
	x[1] += one_minus_t_cubed * p0[1];

	x[0] += 3.0*one_minus_t_squared*t * p1[0];
	x[1] += 3.0*one_minus_t_squared*t * p1[1];

	x[0] += 3.0*one_minus_t*t_squared * p2[0];
	x[1] += 3.0*one_minus_t*t_squared * p2[1];

	x[0] += t_cubed * p3[0];
	x[1] += t_cubed * p3[1];

}

