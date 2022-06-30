
#include <math.h>
#include "vectorops.h" 

float norm(float* n) {	return sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );}

void normalize(float* n, float* r) {
	float nn = norm( n );
	r[0] = n[0] / nn ;
	r[1] = n[1] / nn ;
	r[2] = n[2] / nn ;
}

void cross_product(float* a, float* b, float* n) {
	n[0] = a[1]*b[2] - a[2]*b[1];
	n[1] = a[2]*b[0] - a[0]*b[2];
	n[2] = a[0]*b[1] - a[1]*b[0];
}

void subtract(float* r, float* a, float* b) {
	r[0] = a[0] - b[0] ;
	r[1] = a[1] - b[1] ;
	r[2] = a[2] - b[2] ;
}

void add(float* r, float* a, float* b) {
	r[0] = a[0] + b[0] ;
	r[1] = a[1] + b[1] ;
	r[2] = a[2] + b[2] ;
}

void subtract2(float* r, float* a, float* b) {
	r[0] = a[0] - b[0] ;
	r[1] = a[1] - b[1] ;
}

void add2(float* r, float* a, float* b) {
	r[0] = a[0] + b[0] ;
	r[1] = a[1] + b[1] ;
}

void set2(float* r, float* a) {
	r[0] = a[0] ;
	r[1] = a[1] ;
}

float norm22(float* n) {	return ( n[0]*n[0]+n[1]*n[1] );}

float dist22(float* a, float* b) {	
    return ( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) );
}

float dist2(float* a, float* b) {	
    return sqrt( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) );
}


float weighted_sum2(float* r, float w1, float w2, float* v1, float* v2) {	
    r[0] = w1*v1[0] + w2*v2[0];
    r[1] = w1*v1[1] + w2*v2[1];
}
