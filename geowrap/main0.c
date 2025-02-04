
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <GLFW/glfw3.h>
#include "nanosvg.h"
#include "tess.h"
#include "tesselator.h"
#include "triangles.h"
#include "surface.h"
#include "vectorops.h"

//#define TEST_GRID 0
//#define CONE_ONLY 0
//#define SHELL 0
#define FLAT 0
#define REVERSE_MESH_CONTOURS 0
//#define REVERSE_PATTERN_CONTOURS 1

struct ConeParams {
	float w;
	float HH;
	float foot;
	float mouth;
	float height;
	float max_z;
	float min_z;
	float thickness;
	float funnel_height;
	//-v w=$maxx -v HH=$maxy -v foot=$foot -v mouth=$mouth -v height=$height 
} cone ;


//https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
#define M(I,J) m[I*3+J]
/*
    [ cos(theta)+ux^2*(1-cos(theta))        ux*uy*(1-cos(theta)) - uz*sin(theta)    ux*uz*(1-cos(theta)) + uy*sin(theta)  ]
R = [ uy*ux*(1-cos(theta)) + uz*sin(theta)  cos(theta)+uy^2*(1-cos(theta))          uy*uz*(1-cos(theta)) - ux*sin(theta) ]
    [ uz*ux*(1-cos(theta)) - uy*sin(theta)  uz*uy*(1-cos(theta)) + ux*sin(theta)    cos(theta)+uz^2*(1-cos(theta))        ]
*/

pottoy_spec_t spec;  
void transf_X(float* xyz, float* uv ) {
	faceted_X(xyz,uv, &spec ) ;
};

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

void* stdAlloc(void* userData, unsigned int size)
{
	int* allocated = ( int*)userData;
	TESS_NOTUSED(userData);
	*allocated += (int)size;
	return malloc(size);
}

void stdFree(void* userData, void* ptr)
{
	TESS_NOTUSED(userData);
	free(ptr);
}

struct MemPool
{
	unsigned char* buf;
	unsigned int cap;
	unsigned int size;
};

void* poolAlloc( void* userData, unsigned int size )
{
	//printf("call poolAlloc with %d %d\n",size,(size+0x7) & ~0x7 );

	struct MemPool* pool = (struct MemPool*)userData;
	size = (size+0x7) & ~0x7;
	if (pool->size + size < pool->cap)
	{
		unsigned char* ptr = pool->buf + pool->size;
		pool->size += size;
		//printf("grew pool to %d\n",pool->size);
		return ptr;
	}
	//printf("out of mem: %d < %d!\n", pool->size + size, pool->cap);
	return 0;
}

void poolFree( void* userData, void* ptr )
{
	// empty
	TESS_NOTUSED(userData);
	TESS_NOTUSED(ptr);
}


// Undefine this to see non-interactive heap allocator version.
#define USE_POOL 1


int run = 1;
int cdt = 0;

//a×b=(a2b3−a3b2)i−(a1b3−a3b1)j+(a1b2−a2b1)k.

void idTransform(float* xprime, float *x) {
	xprime[0] = x[0]  ;
	xprime[1] = x[1]  ;
	xprime[2] = x[2]  ;
}

//}


void tangent(float* t, float* rs, float* e, void(*trnsfrm)(float*, float*)) {

   	float rs1[3];
	float rs2[3];

   	float x1[3];
	float x2[3];
	float t0[3];
 
    add(     &rs1[0],rs,e);
    subtract(&rs2[0],rs,e);

    trnsfrm(&x1[0],&rs1[0]);
    trnsfrm(&x2[0],&rs2[0]);

    subtract(&t0[0],&x2[0],&x1[0]);
   // printf("tgt %f %f %f\n",t0[0],t0[1],t0[2]);    fflush(stdout);

    normalize(t0,t);
    //printf("tgt %f %f %f\n",t[0],t[1],t[2]);    fflush(stdout);

}

void surface_norm(float* n, float* rs, void(*trnsfrm)(float*, float*)) {

	float t1[3];
	float t2[3];

    float e1[3] = {0.001,0.0,0.0};
    float e2[3] = {0.0,0.001,0.0};

    tangent(&t1[0],rs,&e1[0],trnsfrm);
    tangent(&t2[0],rs,&e2[0],trnsfrm);

    //printf("t1 %f %f %f\n",t1[0],t1[1],t1[2]);    fflush(stdout);
   // printf("t2 %f %f %f\n",t2[0],t2[1],t2[2]);    fflush(stdout);

	float nn[3];

    cross_product(t1,t2,&nn[0]);
    normalize(nn,n);
}
// take a function transfm (r,s)->(xp,yp,zp), and use it to map (r,s,z)->transf(r,s) + z*normal_at(r,s)
void transform2surf(float* xprime, float *x, void(*trnsfrm)(float*, float*)) {

    float rs[3];
    float p0[3];
    float n[3];
    rs[0] = x[0];
    rs[1] = x[1];
    rs[2] = 0.0;
    trnsfrm(&p0[0],&rs[0]);

    surface_norm(&n[0], &rs[0], trnsfrm );
   // printf("surface norm %f %f %f\n",n[0],n[1],n[2]);    fflush(stdout);
    n[0]*=x[2];
    n[1]*=x[2];
    n[2]*=x[2];

    add(xprime, &p0[0], &n[0]);

}

void transform(float* xprime, float *x) {
	//transf_X
    transform2surf( xprime, x, transf_X );
	//xprime[0] = ( x[0] / cone.w ) * cone.foot;
	//xprime[2] = ( x[1] / cone.w ) * cone.foot;
	//xprime[1] = ( x[2] / cone.w ) * cone.foot;
}

void triangle_normal(float* v1, float* v2, float* v3, float *n) {

	float a[3];
	float b[3];

	a[0] = v2[0] - v1[0];
	a[1] = v2[1] - v1[1];
	a[2] = v2[2] - v1[2];

	b[0] = v3[0] - v1[0];
	b[1] = v3[1] - v1[1];
	b[2] = v3[2] - v1[2];

	cross_product( &a[0], &b[0], &n[0] );

	float l = sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );
	if (l==0.0) { 
		n[0]=0.0;
		n[1]=0.0;
		n[2]=0.0;
		return;
	}
	n[0] /= l;
	n[1] /= l;
	n[2] /= l;
}


void print_triangle_raw( float* v1, float* v2, float* v3 , FILE* fp) {


	float n[3];
	triangle_normal(v1,v2,v3,&n[0]);

	//printf("  facet normal %f %f %f\n", n[0], n[1], n[2]);
	//fprintf(fp,"  facet normal %f %f %f\n", 0.0,0.0,0.0);
    fprintf(fp,"  facet normal %f %f %f\n", n[0], n[1], n[2]);
	fprintf(fp,"    outer loop\n");
	fprintf(fp,"      vertex %f %f %f\n",v1[0],v1[1],v1[2]);
	fprintf(fp,"      vertex %f %f %f\n",v2[0],v2[1],v2[2]);
	fprintf(fp,"      vertex %f %f %f\n",v3[0],v3[1],v3[2]);
	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");


}

void print_triangle( float* z1, float* z2, float* z3 , FILE* fp) {

	float v1[3];
	float v2[3];
	float v3[3];

#if FLAT
	idTransform(v1,z1);
	idTransform(v2,z2);
	idTransform(v3,z3);
#else
	transform(v1,z1);
	transform(v2,z2);
	transform(v3,z3);
#endif

	float n[3];
	triangle_normal(v1,v2,v3,&n[0]);

	if (n[0]==0.0 && n[1]==0.0 && n[2]==0.0) {return;}

	fprintf(fp,"  facet normal %f %f %f\n", n[0], n[1], n[2]);
	fprintf(fp,"    outer loop\n");
	fprintf(fp,"      vertex %f %f %f\n",v1[0],v1[1],v1[2]);
	fprintf(fp,"      vertex %f %f %f\n",v2[0],v2[1],v2[2]);
	fprintf(fp,"      vertex %f %f %f\n",v3[0],v3[1],v3[2]);
	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");

}

void print_slab_triangles(float thickness, FILE* fp) {

	float p1[3] = { 0, 0, -thickness };
	float p2[3] = { 0, cone.HH, -thickness };
	float p3[3] = { cone.w, cone.HH, -thickness };
	float p4[3] = { cone.w, 0, -thickness };

	print_triangle(&p1[0],&p3[0],&p2[0], fp);
	print_triangle(&p1[0],&p4[0],&p3[0], fp);

	float p1b[3] = { 0, 0, -2.0*thickness };
	float p2b[3] = { 0, cone.HH, -2.0*thickness };
	float p3b[3] = { cone.w, cone.HH, -2.0*thickness };
	float p4b[3] = { cone.w, 0, -2.0*thickness };

	print_triangle(&p1b[0],&p2b[0],&p3b[0], fp);
	print_triangle(&p1b[0],&p3b[0],&p4b[0], fp);

	print_triangle(&p1b[0],&p4[0],&p1[0], fp);
	print_triangle(&p1b[0],&p4b[0],&p4[0], fp);

	print_triangle(&p4b[0],&p3[0],&p4[0], fp);
	print_triangle(&p4b[0],&p3b[0],&p3[0], fp);

	print_triangle(&p2[0],&p3[0],&p2b[0], fp);
	print_triangle(&p2b[0],&p3[0],&p3b[0], fp);

	print_triangle(&p2b[0],&p1[0],&p2[0], fp);
	print_triangle(&p2b[0],&p1b[0],&p1[0], fp);

}

float dot(float* x, float* y) {
	return ( x[0]*y[0]+x[1]*y[1]+x[2]*y[2] );
}

void print_patch( float* x10,float* x20,float* x30, float thickness, float chamfer, FILE* fp) {

	float x1[3];// = { x10[0], x10[1], 0 };
	float x2[3];// = { x20[0], x20[1], 0 };
	float x3[3];// = { x20[0], x20[1], 0 };
    float xa[3];
    float xb[3];

	if (vequal2(x10,x20)) {printf("x10==x20 nan\n");}
	if (vequal2(x20,x30)) {printf("x20==x30 nan\n");}

    transform(&x1[0],x10);
    transform(&x2[0],x20);
    transform(&x3[0],x30);

	if (vequal3(x1,x2)) { return;}
	if (vequal3(x2,x3)) { return;}

	if (isnan(x1[0])) {printf("x1[0] is nan\n");}
	if (isnan(x2[0])) {printf("x2[0] is nan\n");}
	if (isnan(x3[0])) {printf("x3[0] is nan\n");}

	float t1[3] = { x2[0] - x1[0], x2[1] - x1[1] , x2[2] - x1[2] };
	float t2[3] = { x3[0] - x2[0], x3[1] - x2[1] , x3[2] - x2[2] };

   // xa[0]=0.5*( x10[0]+x20[0] );
   // xa[1]=0.5*( x10[1]+x20[1] );
   // xa[2]=0.5*( x10[2]+x20[2] );

   // xb[0]=0.5*( x20[0]+x30[0] );
   // xb[1]=0.5*( x20[1]+x30[1] );
  //  xb[2]=0.5*( x20[2]+x30[2] );
    //printf("xa: %f %f %f \n",xa[0],xa[1],xa[2]); //,x2[0],x2[1],x2[2]);

	float normal2[3]; // = { 0,0, -1 };
//	float normalB[3]; // = { 0,0, -1 };

    surface_norm(&normal2[0], &x20[0], transf_X );
//    surface_norm(&normalB[0], &xb[0], X );
	if (isnan(normal2[0])) {printf("normal2[0] is nan\n");}
	//if (isnan(x2[0])) {printf("x2[0] is nan\n");}
	//if (isnan(x3[0])) {printf("x3[0] is nan\n");}

//	float normal[3] = { 0,0, -1 };

	float tn1[3];
	float tn2[3];

	float m1[9];
	float m2[9];

	//float t1xt2[3];

	normalize(&t1[0],&tn1[0]);
	normalize(&t2[0],&tn2[0]);
	
	if (isnan(tn1[0])) {printf("tn1[0] is nan\n");}
	if (isnan(tn2[0])) {printf("tn2[0] is nan\n");}


	//cross_product(&tn1[0],&tn2[0],&t1xt2[0]);

	//if (dot(&t1xt2[0], normal ) <0.0) {return ; }

	rotationMatrix(chamfer, tn1, &m1[0]);
	rotationMatrix(chamfer, tn2, &m2[0]);

	float chamfer_v1[3]; //  = { 0,0, 1 };
	float chamfer_v2[3]; //  = { 0,0, 1 };

	matrixMultiply(&m1[0],&normal2[0],&chamfer_v1[0]);
	matrixMultiply(&m2[0],&normal2[0],&chamfer_v2[0]);

	if (isnan(chamfer_v1[0])) {printf("chamfer_v1[0] is nan\n");}
	if (isnan(chamfer_v2[0])) {printf("chamfer_v2[0] is nan\n");}

	float r = thickness / cos(chamfer);
    //printf("r: %f \t %f \t %f\n",r,thickness,chamfer);
	float v1[3] = { x2[0], x2[1], x2[2] } ;
	//float v2[3] = { this_x, this_y, 0 } ;
	float v2[3] = { x2[0] - r*chamfer_v1[0], x2[1] - r*chamfer_v1[1], x2[2] - r*chamfer_v1[2]  } ;
	float v3[3] = { x2[0] - r*chamfer_v2[0], x2[1] - r*chamfer_v2[1], x2[2] - r*chamfer_v2[2]  } ;
	float v4[3] = { x2[0] - thickness*normal2[0], x2[1] - thickness*normal2[1], x2[2] - thickness*normal2[2] } ;

	print_triangle_raw(v1,v3,v2, fp); 
	print_triangle_raw(v2,v3,v4, fp); 

}

void print_wedge(float* x10, float* x20,float thickness, float chamfer, FILE* fp ) {

    //printf("x10: %f %f %f ; x20: %f %f %f \n",x10[0],x10[1],x10[2],x20[0],x20[1],x20[2]);
    
	float x1[3];// = { x10[0], x10[1], 0 };
	float x2[3];// = { x20[0], x20[1], 0 };
    float xa[3];

    transform(&x1[0],x10);
    transform(&x2[0],x20);

	if (vequal3(x1,x2)) { return;}
	//if (vequal3(x2,x3)) {printf("x2==x3 nan\n"); return;}

   // printf("x1: %f %f %f ; x2: %f %f %f \n",x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);

	float t[3] = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] };  // tangent
    /*
    xa[0]=0.5*( x10[0]+x20[0] );
    xa[1]=0.5*( x10[1]+x20[1] );
    xa[2]=0.5*( x10[2]+x20[2] );
    printf("xa: %f %f %f \n",xa[0],xa[1],xa[2]); //,x2[0],x2[1],x2[2]);
    */
	//float normal[3]; // = { 0,0, -1 };
	float normal1[3]; // = { 0,0, -1 };
	float normal2[3]; // = { 0,0, -1 };

    surface_norm(&normal1[0], &x10[0], transf_X );
    surface_norm(&normal2[0], &x20[0], transf_X );

//    surface_norm(&normal[0], &xa[0], X );
 //   printf("norm: %f %f %f \n",normal[0],normal[1],normal[2]); //,x2[0],x2[1],x2[2]);

	float chamfer_v1[3] ;// = { normal[0], normal[1], normal[2] };
	float chamfer_v2[3] ;// = { normal[0], normal[1], normal[2] };
	float tn[3] ;
	float m[9];
	normalize(&t[0],&tn[0]);  
	rotationMatrix(chamfer, tn, &m[0]);
	matrixMultiply(&m[0],&normal1[0],&chamfer_v1[0]);
	matrixMultiply(&m[0],&normal2[0],&chamfer_v2[0]);
   // printf("chamfer_v: %f %f %f \n",chamfer_v[0],chamfer_v[1],chamfer_v[2]); //,x2[0],x2[1],x2[2]);

	float r = thickness / cos(chamfer);
    //r = thickness;
	//float v1[3] = { last_x, last_y, 0 } ;
	//float v2[3] = { this_x, this_y, 0 } ;

	float x3[3] = { x1[0] - r*chamfer_v1[0], x1[1] - r*chamfer_v1[1], x1[2] - r*chamfer_v1[2]  } ;
	float x4[3] = { x2[0] - r*chamfer_v2[0], x2[1] - r*chamfer_v2[1], x2[2] - r*chamfer_v2[2]  } ;

	float x5[3] = { x1[0] - thickness*normal1[0] , x1[1] - thickness*normal1[1] ,x1[2]  - thickness*normal1[2] };
	float x6[3] = { x2[0] - thickness*normal2[0] , x2[1] - thickness*normal2[1] ,x2[2]  - thickness*normal2[2] };

//	v1 = {}; 

	print_triangle_raw(x1,x2,x3, fp); 
	print_triangle_raw(x3,x2,x4, fp); 

	//print_triangle(x1,x3,x5); 
	//print_triangle(x4,x2,x6); 
    
	print_triangle_raw(x6,x5,x3, fp); 
	print_triangle_raw(x4,x6,x3, fp); 
	print_triangle_raw(x2,x1,x5, fp); 
	print_triangle_raw(x2,x5,x6, fp); 
    

}

void print_quad(float last_x, float last_y, float this_x, float this_y, float thickness, FILE* fp) {

	float v1[3] = { last_x, last_y, 0.0 } ;
	float v2[3] = { this_x, this_y, 0.0 } ;
	float v3[3] = { last_x , last_y , -thickness  } ;
	float v4[3] = { this_x , this_y , -thickness  } ;


	print_triangle(v1,v3,v2, fp); 
	print_triangle(v3,v4,v2, fp); 

}


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

int main(int argc, char *argv[])
{
	//struct SVGPath* bg;
	NSVGimage *bg = NULL;

	//struct SVGPath* it;
	float bounds[4];
	float chamfer;
	//coneParams cone;

	int width,height,i,j;
	//float fb[8]  ;

	float cx,cy;
	float expansion_fact = 1.0;
//	float expansion_fact = argc >= 11 ? 1.0/(1.0-atof(argv[10])) : 1.0 ;
//	printf("expansion_fact  %f\n",expansion_fact); fflush(stdout);

	//float t = 0.0f, pt = 0.0f;
	TESSalloc ma;
	TESStesselator* tess = 0;
	const int nvp = 3;
	unsigned char* vflags = 0;

	struct MemPool pool;
	unsigned char* mem;   // [1024*1024*20];
	int nvflags = 0;

    if (argc==2) {
        write_surface_obj(argv[1]);     
        //FILE* stlfile = fopen(argv[3], "wt"); 
        exit(0);
    }
		//printf("argc=%d\n",argc);

	int n_sectors = 90;
	int n_levels = 50;

    if (argc==3) {
		read_spec(argv[1], &spec );
	} else {
		read_spec(argv[4], &spec );
	}

	if (spec.facet) {
		int spf = n_sectors / spec.n_facets ;
		if (spf < 10 ) {spf = 10;}
		n_sectors = spec.n_facets * (spf + 1);
		printf("n_sectors = %d\n", n_sectors );
	}
		
    if (argc==3) {
		printf("read spec & write obj.\n");
		read_spec(argv[1], &spec );

		float bbox[6];
		surface_obj_bbox(&bbox[0], &spec) ;
		printf("bbox: x:  %f - %f\n", bbox[0], bbox[3]);
		printf("bbox: y:  %f - %f\n", bbox[1], bbox[4]);
		printf("bbox: z:  %f - %f\n", bbox[2], bbox[5]);


		//int n_sectors = 90;


  		write_surface_obj2(argv[2], &spec, n_sectors, n_levels);     



    	//write_surface_obj2(argv[2], &spec);     

        exit(0);
    }

    if (argc<6) {
        printf("usage: wemboss <input.svg> <input.obj> <out.stl> <surface.json> <thickness> <chamfer>\n");
        return(0);
    }

    FILE* stlfile = fopen(argv[3], "wt"); 
	fprintf(stlfile,"solid x\n");


	//printf("hi %s %s\n",argv[0],argv[1]); fflush(stdout);
	//bg = svgParseFromFile(argv[1]);
	bg = nsvgParseFromFile(argv[1], "px", 98.0f);
	if (!bg) { printf("error parsing %s\n",argv[1]); return -1; }
	
    //parse .obj file to get the interpolation mesh
    //bg = nsvgParseFromFile(argv[2], "px", 98.0f);
    printf("parse file %s\n",argv[2]);

	float thickness = atof(argv[5])*expansion_fact;
	cone.thickness = thickness;
	chamfer = (M_PI*atof(argv[6])/180.0);
	float scale=1.0;
	if (argc>=8) { 
		scale = atof(argv[7]); 
		scale_spec(scale,&spec);
	}
	float shell_thickness = 1.0;
	if (argc>=9) { 
		shell_thickness = atof(argv[8]); 
	}

#ifdef TEST_GRID
	float chamfer0 = chamfer;
#endif

	int np = 0;
	int done_init=0;

	float *p0, *p1, *p2, *p3; 

	int n_points = 0;
	int n_paths = 0;
    int ppi = 0;

#define SAMPLES_PER_BEZ 5

	for (NSVGshape *shape = bg->shapes; shape != NULL; shape = shape->next) {
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {	
			n_paths += 1;	
			for (int i = 0; i < path->npts-1; i+=3) {
				n_points += SAMPLES_PER_BEZ;
                for (float t=0.0; t<1.0; t+=1.0/((float) (SAMPLES_PER_BEZ+1) ) ) {
                    ppi+=2;
                }
			}
		}
	}

	int *path_lengths = malloc( n_paths * sizeof(int) );
	int *path_offsets = malloc( n_paths * sizeof(int) );
	float *path_points = malloc( n_points * 2 * sizeof(float) );

    for (int i=0; i<n_points; i++) {
        path_points[2*i]  =0; 
        path_points[2*i+1]=0; 
    }


	int pi = 0;
	 ppi = 0;

	// read in the curves from the SVG and compute bounding box.
	for (NSVGshape *shape = bg->shapes; shape != NULL; shape = shape->next) {
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {
			//printf("x:\nx:\n");
			path_offsets[pi]=ppi;
			path_lengths[pi]=0;
			for (int i = 0; i < path->npts-1; i+=3) {

				//float* p = &path->pts[i*2];

				p0 = &path->pts[i*2];
				p1 = &path->pts[((i+1)%path->npts)*2];
				p2 = &path->pts[((i+2)%path->npts)*2];
				p3 = &path->pts[((i+3)%path->npts)*2];

				float xx[2] = {0,0};
                float dt = 1.0 / ( (float) SAMPLES_PER_BEZ );
                float t=-dt;
				for ( int ti=0; ti<SAMPLES_PER_BEZ ; ti++ ) {
                    t+=dt;
					path_lengths[pi]+=1;
					cubicBez(t,p0,p1,p2,p3,&xx[0]);

					const float x = xx[0];
					const float y = xx[1];

					if (!done_init) {
						done_init = 1;
						//printf("init\n");
							bounds[0] = bounds[2] = x;
							bounds[1] = bounds[3] = y;
					} 

					np+=1;
					if (x < bounds[0]) bounds[0] = x;
					if (y < bounds[1]) bounds[1] = y;
					if (x > bounds[2]) bounds[2] = x;
					if (y > bounds[3]) bounds[3] = y;
					//printf( "x: %f %f\n", p[0],p[1] );

					//printf( "x: %f %f\n", xx[0],xx[1] );
					path_points[ppi  ] = xx[0];
					path_points[ppi+1] = xx[1];
					ppi+=2;

				}

				//drawCubicBez(p[0],p[1], p[2],p[3], p[4],p[5], p[6],p[7]);
			}
			pi+=1;
		}
	
	}
    if(ppi>n_points*2) {printf("%d >= %d\n",ppi,n_points*2); exit(1);}
    //return(0);

	cx = (bounds[0]+bounds[2])/2;
	cy = (bounds[1]+bounds[3])/2;
	width = bounds[2]-bounds[0]; 
	height = bounds[3]-bounds[1];

	cone.w = width;
	cone.HH = bounds[3]-bounds[1]-1.0;

    printf("width: %f\n",(float) width);
    MeshTriangles* mt = parse_triangles(argv[2],(float) width);
	fflush(stdout);

    for (int i=0; i<n_points; i++) {
       // printf("interpolate %d %d\n",i,n_points  );
        mesh_interpolation(mt, & path_points[2*i], & path_points[2*i] );
    }

    printf("did mesh interpolation\n");
	fflush(stdout);
	np += n_sectors * n_levels * 9 * 3;
	mem = malloc( np*4096 );
	//printf("mem=%x\n",mem);
	memset(mem,0, np*4096 );
	printf("allocate %d K\n",np*4096/1024); fflush(stdout);

	pool.size = 0;
	pool.cap = np*4096; //sizeof(mem);
	pool.buf = mem;
	memset(&ma, 0, sizeof(ma));
	//memset(&ma,0, sizeof(ma))
	ma.memalloc = poolAlloc;
	ma.memfree = poolFree;
	ma.userData = (void*)&pool;
	ma.extraVertices = 4096; // realloc not provided, allow 256 extra vertices.


		pool.size = 0; // reset pool
		tess = tessNewTess(&ma);
		if (tess)
		{
			tessSetOption(tess, TESS_CONSTRAINED_DELAUNAY_TRIANGULATION, 0);

//	for (NSVGshape *shape = bg->shapes; shape != NULL; shape = shape->next) {
//		for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {
//			printf("x:\nx:\n");
//			for (int i = 0; i < path->npts; ++i) {
	
#if REVERSE_PATTERN_CONTOURS
                tessSetOption(tess,TESS_REVERSE_CONTOURS,1);
#endif

    	for (int i =0; i<n_paths; i++) {
			tessAddContour(tess, 2, &path_points[ path_offsets[i] ] , sizeof(float)*2, path_lengths[i]);
		}
#if REVERSE_PATTERN_CONTOURS
                tessSetOption(tess,TESS_REVERSE_CONTOURS,1);
#endif



/*
                // add the mesh trianlges as paths
                for (int i=0; i<mt->ntriangles; i++) {
                    tessAddContour(tess, 2, (void*) &(mt->paths[i*6]), sizeof(float)*2, 3);

                    
                }
                */

		/*
			for (NSVGshape *shape = bg->shapes; shape != NULL; shape = shape->next) {
				for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {
					tessAddContour(tess, 2, path->pts, sizeof(float)*2, path->npts);
				}
			}
			*/


//				for (int i = 0; i < path->npts-1; i += 3) {

			//for (it = bg; it != NULL; it = it->next)
//				tessAddContour(tess, 2, it->pts, sizeof(float)*2, it->npts);
			
		
#if 1
		// First combine contours and then triangulate, this removes unnecessary inner vertices.
			if (tessTesselate(tess, TESS_WINDING_ODD, TESS_BOUNDARY_CONTOURS, 0, 0, 0))
			//if (tessTesselate(tess, TESS_WINDING_ABS_GEQ_TWO, TESS_BOUNDARY_CONTOURS, 0, 0, 0))
			{
				const float* verts = tessGetVertices(tess);
				const int* vinds = tessGetVertexIndices(tess);
				const int nverts = tessGetVertexCount(tess);
				const int* elems = tessGetElements(tess);
				const int nelems = tessGetElementCount(tess);

                printf("n elements: %d\n",nelems);

				if (nverts > nvflags)
				{
					if (vflags)
						free(vflags);
					nvflags = nverts;
				//	printf("alloc %d vflags",nvflags);

					vflags = (unsigned char*)malloc(sizeof(unsigned char)*nvflags);
				}

				if (vflags)
				{
					// Vertex indices describe the order the indices were added and can be used
					// to map the tesselator output to input. Vertices marked as TESS_UNDEF
					// are the ones that were created at the intersection of segments.
					// That is, if vflags is set it means that the vertex comes from intersegment.
					for (i = 0; i < nverts; ++i)
						vflags[i] = vinds[i] == TESS_UNDEF ? 1 : 0;
				}


                printf("nelems %d\n",nelems);  fflush(stdout);
				for (i = 0; i < nelems; ++i)
				{
					int b = elems[i*2];
					int n = elems[i*2+1];

                    
						//printf("tc:\n");
						//printf("tc:\n");
				//	printf("add contour %d/%d %d %d\n",i,nelems-1,b,n);


					//for (j = 0; j < n; ++j)
					//{
					//	printf("tc: %f %f\n",verts[b*2+j*2], verts[b*2+j*2+1]);
				//	}

					
		

					tessAddContour(tess, 2, &verts[b*2], sizeof(float)*2, n);
				}


#if REVERSE_MESH_CONTOURS
                tessSetOption(tess,TESS_REVERSE_CONTOURS,1);
#endif
//                tessSetOption(tess,TESS_REVERSE_CONTOURS,1);

				add_triangle_contours(tess, &spec, 3*n_sectors, 3*n_levels) ;
				
				/*
                // add the mesh trianlges as paths
                for (int i=0; i<mt->npaths; i++) {
                //    tessAddContour(tess, 2, (void*) &(mt->paths[i*6]), sizeof(float)*2, 3);
					printf("# %f,%f %f,%f %f,%f \n",
						mt->paths[i*6],
						mt->paths[i*6+1],
						mt->paths[i*6+2],
						mt->paths[i*6+3],
						mt->paths[i*6+4],
						mt->paths[i*6+5]
					);

                }*/
				
#if REVERSE_MESH_CONTOURS
				tessSetOption(tess,TESS_REVERSE_CONTOURS,0);
#endif
                
                // add rectangular strips in image coordinates.  
				// These form the object mesh...  

    /*
				float dx = width / 200.0;
				for (float x=bounds[0]-1.0;x<bounds[2];x+=dx) {
					fb[0] = x;
					fb[1] = bounds[1]-1.0;
					fb[2] = x;
					fb[3] = bounds[3]+1.0;
					fb[4] = x+dx;
					fb[5] = bounds[3]+1.0;
					fb[6] = x+dx;
					fb[7] = bounds[1]-1.0;
					tessAddContour(tess, 2, fb, sizeof(float)*2, 4);

                    printf("tt:\ntt:\n");
                    for (j = 0; j < 4; ++j)
					{
						printf("tt: %f %f\n",fb[j*2], fb[j*2+1]);
					}

				}
                */

			} else { printf("tessTesslate returned zero A\n"); exit(0); }


			if (tessTesselate(tess, TESS_WINDING_ABS_GEQ_TWO, TESS_BOUNDARY_CONTOURS, 0, 0, 0)) {

				const float* verts = tessGetVertices(tess);
				const int* vinds = tessGetVertexIndices(tess);
				const int nverts = tessGetVertexCount(tess);
				const int* elems = tessGetElements(tess);
				const int nelems = tessGetElementCount(tess);
                printf("n elements: %d\n",nelems);

				if (nverts > nvflags)
				{
					if (vflags)
						free(vflags);
					nvflags = nverts;
				//	printf("alloc %d vflags",nvflags);

					vflags = (unsigned char*)malloc(sizeof(unsigned char)*nvflags);
				}

				if (vflags)
				{
					// Vertex indices describe the order the indices were added and can be used
					// to map the tesselator output to input. Vertices marked as TESS_UNDEF
					// are the ones that were created at the intersection of segments.
					// That is, if vflags is set it means that the vertex comes from intersegment.
					for (i = 0; i < nverts; ++i)
						vflags[i] = vinds[i] == TESS_UNDEF ? 1 : 0;
				}
				for (i = 0; i < nelems; ++i)
				{
					int b = elems[i*2];
					int n = elems[i*2+1];
					//	printf("tc:\n");
					//	printf("tc:\n");
					//printf("add contour %d/%d %d %d\n",i,nelems-1,b,n);

					/*for (j = 0; j < n; ++j)
					{
						printf("tc: %f %f\n",verts[b*2+j*2], verts[b*2+j*2+1]);
					}*/

					tessAddContour(tess, 2, &verts[b*2], sizeof(float)*2, n);
				}

				//side walls 
				for (i = 0; i < nelems; ++i)
				{
					int b = elems[i*2];
					int n = elems[i*2+1];

					//	printf("tc:\n");
					//	printf("tc:\n");

				//	float last_x = verts[b*2];
				//	float last_y = verts[b*2+1];

					//float this_x; float this_y;
					for (j = 0; j < n; ++j)
					{
						float x1[3] = { verts[b*2+(j%n)*2], verts[b*2+(j%n)*2+1],0 } ;
						float x2[3] = { verts[b*2+((j+1)%n)*2], verts[b*2+((j+1)%n)*2+1],0 } ;
					
						//this_x = verts[b*2+j*2];
						//this_y = verts[b*2+j*2+1];

#ifdef TEST_GRID 
					    float theta1 = 2.0*M_PI*( x1[0] / cone.w );
						int sector1 = ((int) (theta1 / (2.0*M_PI/2.0)))%2; 
						chamfer = chamfer0 * ((float) sector1);
						//printf("chamfer %d %f %f\n",sector1,chamfer,0.0);//180.0*(theta1/M_PI));
#endif

#ifndef CONE_ONLY
						if ( dist2(&x2[0],&x1[0])>0 ) {
#if 1
						
						if (fabs(chamfer)>0.001) {
						//	printf("chamfer %f\n",chamfer) ;
							print_wedge(&x2[0],&x1[0],1.1*thickness,chamfer,stlfile);
						} 
#endif					
						print_quad( x1[0], x1[1], x2[0], x2[1], 1.1*thickness, stlfile);
						}
#endif
					}
				}


				//side wall patches between chamfered rectangles
				for (i = 0; i < nelems; ++i)
				{
				
					int b = elems[i*2];
					int n = elems[i*2+1];

					for (j = 0; j < n; ++j)
					{
						float x1[3] = { verts[b*2+(j%n)*2], verts[b*2+(j%n)*2+1] ,0} ;
						float x2[3] = { verts[b*2+((j+1)%n)*2], verts[b*2+((j+1)%n)*2+1] ,0} ;
						float x3[3] = { verts[b*2+((j+2)%n)*2], verts[b*2+((j+2)%n)*2+1] ,0} ;
						//float x4[2] = { verts[b*2+((j+3)%n)*2], verts[b*2+((j+3)%n)*2+1] } ;

#ifdef NAN
				if (x1[0]==NAN) {printf("x1[0]=NAN\n");}
				if (x1[1]==NAN) {printf("x1[0]=NAN\n");}
				if (x1[2]==NAN) {printf("x1[0]=NAN\n");}
				if (x2[0]==NAN) {printf("x1[0]=NAN\n");}
				if (x2[1]==NAN) {printf("x1[0]=NAN\n");}
				if (x2[2]==NAN) {printf("x1[0]=NAN\n");}
				if (x3[0]==NAN) {printf("x1[0]=NAN\n");}
				if (x3[1]==NAN) {printf("x1[0]=NAN\n");}
				if (x3[2]==NAN) {printf("x1[0]=NAN\n");}
#endif

#ifdef TEST_GRID 
					    float theta1 = 2.0*M_PI*( x2[0] / cone.w );
						int sector1 = ((int) (theta1 / (2.0*M_PI/2.0)))%2; 
						chamfer = chamfer0 * ((float) sector1);
#endif
#ifndef CONE_ONLY

						if (fabs(chamfer)>0.001) {
							//if (dist22(&x2[0],&x1[0])==0.0) {fprintf(stlfile,"x2==x1\n");}
							//if (dist22(&x2[0],&x3[0])==0.0) {fprintf(stlfile,"x2==x3\n");}
							//fprintf(stlfile,"patch %d %d\n",i,j);
							print_patch( x3,x2,x1, 1.1*thickness, chamfer, stlfile);
						}
#endif
					}

				}


			//for (it = bg; it != NULL; it = it->next)
			//	tessAddContour(tess, 2, it->pts, sizeof(float)*2, it->npts);
			
            /*
				float dx = width / 200.0;
				for (float x=bounds[0]-1.0;x<bounds[2];x+=dx) {
					fb[0] = x;
					fb[1] = bounds[1]-1.0;
					fb[2] = x;
					fb[3] = bounds[3]+1.0;
					fb[4] = x+dx;
					fb[5] = bounds[3]+1.0;
					fb[6] = x+dx;
					fb[7] = bounds[1]-1.0;
					tessAddContour(tess, 2, fb, sizeof(float)*2, 4);

				}*/
#if REVERSE_MESH_CONTOURS
				tessSetOption(tess,TESS_REVERSE_CONTOURS,1);
#endif
                // add the mesh trianlges as paths

				add_triangle_contours(tess, &spec, 3*n_sectors, 3* n_levels) ;
/*
                for (int i=0; i<mt->npaths; i++) {
                    tessAddContour(tess, 2, (void*) &(mt->paths[i*6]), sizeof(float)*2, 3);
                  
                }
				*/
#if REVERSE_MESH_CONTOURS
                tessSetOption(tess,TESS_REVERSE_CONTOURS,0);
#endif

				//printf("done adding contours\n"); //TESS_WINDING_POSITIVE,
			if (!tessTesselate(tess, TESS_WINDING_ABS_GEQ_TWO, TESS_POLYGONS, nvp, 2, 0)) {
				//if (!tessTesselate(tess,TESS_WINDING_ABS_GEQ_TWO, TESS_POLYGONS, nvp, 2, 0))
					tess = 0;
					printf("tesselate returned zero 2nd time\n");
					fflush(stdout);
					//printf("tesselate returned zero 2nd time; out of memory: %d\n",tess->outOfMemory);
				}
			}

#else
			if (tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS, nvp, 2, 0)) {
				printf("tesselated");
			}
#endif
			else {
				printf("tesselate returned zero 1st time\n");
				tess = 0;
			}

			//tessTesselate(tess, TESS_WINDING_ABS_GEQ_TWO, TESS_POLYGONS, nvp, 2, 0);
			//tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS, nvp, 2, 0);

			if (tess)
				{
					const float* verts = tessGetVertices(tess);
					//const int* vinds = tessGetVertexIndices(tess);
					const int* elems = tessGetElements(tess);
					//const int nverts = tessGetVertexCount(tess);
					const int nelems = tessGetElementCount(tess);

                    printf("final tesselation n elements: %d\n",nelems);
					// Draw polygons.
					//glColor4ub(255,255,255,128);
					// top
#ifndef CONE_ONLY

//The triangles filling in the mesa tops and bottoms
					for (i = 0; i < nelems; ++i)
					{

						const int* p = &elems[i*nvp];

						float v1[3] = { verts[p[0]*2],verts[p[0]*2+1],0.0*thickness };
						float v2[3] = { verts[p[1]*2],verts[p[1]*2+1],0.0*thickness };
						float v3[3] = { verts[p[2]*2],verts[p[2]*2+1],0.0*thickness };

						print_triangle(v1,v2,v3, stlfile);

						v1[2] = -1.1*thickness;
						v2[2] = -1.1*thickness;
						v3[2] = -1.1*thickness;

						print_triangle(v1,v3,v2, stlfile);

					}
#endif

					// print out the shell of the surface itself -- outside
                    for (int i=0; i<mt->npaths; i++) {
                        //tessAddContour(tess, 2, &path_points[ path_offsets[i] ] , sizeof(float)*2, path_lengths[i]);
                        float v1[3] = { mt->paths[i*6+0] ,mt->paths[i*6+1] ,-1.0*thickness };
						float v2[3] = { mt->paths[i*6+2] ,mt->paths[i*6+3] ,-1.0*thickness };
						float v3[3] = { mt->paths[i*6+4] ,mt->paths[i*6+5] ,-1.0*thickness };
                        print_triangle(v1,v2,v3, stlfile);
                    }
					
					
					// inside
                    for (int i=0; i<mt->npaths; i++) {
                        //tessAddContour(tess, 2, &path_points[ path_offsets[i] ] , sizeof(float)*2, path_lengths[i]);
                        float v1[3] = { mt->paths[i*6+0] ,mt->paths[i*6+1] ,-1.0*thickness - shell_thickness };
						float v2[3] = { mt->paths[i*6+2] ,mt->paths[i*6+3] ,-1.0*thickness - shell_thickness };
						float v3[3] = { mt->paths[i*6+4] ,mt->paths[i*6+5] ,-1.0*thickness - shell_thickness };
                        print_triangle(v1,v3,v2, stlfile);
                    }
					
				}
		}
	
#if 1
#if FLAT
	//print_slab_triangles(thickness);
#else
	//print_cone_triangles(thickness);
#endif	
#endif
	
	fprintf(stlfile, "endsolid x\n");


    fclose(stlfile);
	return 0;
}
