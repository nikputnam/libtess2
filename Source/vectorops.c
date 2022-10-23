
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

void set_triangle(float* r, float* a, float* b, float* c) {
	set( r,    a );
	set( &r[3],b );
	set( &r[6],c );
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

void diff3(float* x, float* y, float* z) {
    z[0] = x[0]-y[0];
    z[1] = x[1]-y[1];
    z[2] = x[2]-y[2];
}

float dot(float* x, float* y) {
	return ( x[0]*y[0]+x[1]*y[1]+x[2]*y[2] );
}

#define M(I,J) m[I*3+J]

void rotationMatrix(float theta, float* u, float* m) {
	float ux = u[0];
	float uy = u[1];
	float uz = u[2];

    float nn = norm(u);
    /*
    if (fabs(nn-1.0)>0.0000001) {
        printf("nn %f\n",nn); fflush(stdout);
        assert(fabs(nn-1.0)<=0.00001);
    }*/
	M(0,0) = cos(theta)+ux*ux*(1-cos(theta));
	M(0,1) = ux*uy*(1-cos(theta)) - uz*sin(theta);
	M(0,2) = ux*uz*(1-cos(theta)) + uy*sin(theta) ;

	M(1,0) = uy*ux*(1-cos(theta)) + uz*sin(theta);
	M(1,1) = cos(theta)+uy*uy*(1-cos(theta));
	M(1,2) = uy*uz*(1-cos(theta)) - ux*sin(theta) ;

	M(2,0) = uz*ux*(1-cos(theta)) - uy*sin(theta);
	M(2,1) = uz*uy*(1-cos(theta)) + ux*sin(theta) ;
	M(2,2) = cos(theta)+uz*uz*(1-cos(theta)) ;



	if (isnan(M(0,0))) {printf("M(0,0) is nan; theta=%f ux=%f uy=%f uz=%f\n",theta,ux,uy,uz);}
	if (isnan(M(0,1))) {printf("M(0,1) is nan\n");}
	if (isnan(M(0,2))) {printf("M(0,2) is nan\n");}
	if (isnan(M(1,0))) {printf("M(1,0) is nan\n");}
	if (isnan(M(1,1))) {printf("M(1,1) is nan\n");}
	if (isnan(M(1,2))) {printf("M(1,2) is nan\n");}
	if (isnan(M(2,0))) {printf("M(2,0) is nan\n");}
	if (isnan(M(2,1))) {printf("M(2,1) is nan\n");}
	if (isnan(M(2,2))) {printf("M(2,2) is nan\n");}

}

void matrixMultiply(float* m, float* x, float* r) {
	r[0] = M(0,0)*x[0] + M(0,1)*x[1] + M(0,2)*x[2];
	r[1] = M(1,0)*x[0] + M(1,1)*x[1] + M(1,2)*x[2];
	r[2] = M(2,0)*x[0] + M(2,1)*x[1] + M(2,2)*x[2];
}
