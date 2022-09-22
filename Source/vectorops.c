
#include <math.h>
#include "vectorops.h" 
#include <stdio.h>
#include <stdlib.h>

float norm(float* n) {	return sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );}

int vequal2(float* a, float* b) { return ( a[0]==b[0] && a[1]==b[1] ); }
int vequal3(float* a, float* b) { return ( a[0]==b[0] && a[1]==b[1] && a[2]==b[2] ); }

void normalize(float* n, float* r) {
	float nn = norm( n );
	if (nn==0.0) {printf("nn=0.0\n");}
	r[0] = n[0] / nn ;
	r[1] = n[1] / nn ;
	r[2] = n[2] / nn ;

	//if (isnan(r[0])) {printf("is nan\n");}

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


void scale(float* p, float a) {
	p[0] *= a ;
	p[1] *= a ;
	p[2] *= a ;
}


void set(float* r, float* a) {
	r[0] = a[0] ;
	r[1] = a[1] ;
	r[2] = a[2] ;
}


float norm22(float* n) {	return ( n[0]*n[0]+n[1]*n[1] );}

float dist22(float* a, float* b) {	
    return ( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) );
}

float dist2(float* a, float* b) {	
    return sqrt( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) );
}

void weighted_sum3(float* r, float w1, float w2, float* v1, float* v2) {	
    r[0] = w1*v1[0] + w2*v2[0];
    r[1] = w1*v1[1] + w2*v2[1];
    r[2] = w1*v1[2] + w2*v2[2];
}

void weighted_sum2(float* r, float w1, float w2, float* v1, float* v2) {	
    r[0] = w1*v1[0] + w2*v2[0];
    r[1] = w1*v1[1] + w2*v2[1];
}

void weighted_sum2_4(float* r, float w1, float w2, float w3, float w4, float* v1, float* v2, float* v3, float* v4) {	
    r[0] = w1*v1[0] + w2*v2[0] + w3*v3[0] + w4*v4[0] ;
    r[1] = w1*v1[1] + w2*v2[1] + w3*v3[1] + w4*v4[1];
}


float dot2(float* x, float* y){
    return (x[0]*y[0] + x[1]*y[1]);
}

void diff(float* x, float* y, float* z) {
    z[0] = x[0]-y[0];
    z[1] = x[1]-y[1];
}

float dot(float* x, float* y) {
	return ( x[0]*y[0]+x[1]*y[1]+x[2]*y[2] );
}