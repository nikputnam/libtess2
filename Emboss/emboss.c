
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GLFW/glfw3.h>
#include "nanosvg.h"
#include "tesselator.h"



struct ConeParams {
	float w;
	float HH;
	float foot;
	float mouth;
	float height;
	float max_z;
	float min_z;
	float thickness;
	//-v w=$maxx -v HH=$maxy -v foot=$foot -v mouth=$mouth -v height=$height 
} cone ;


//https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
#define M(I,J) m[I*3+J]
/*
    [ cos(theta)+ux^2*(1-cos(theta))        ux*uy*(1-cos(theta)) - uz*sin(theta)    ux*uz*(1-cos(theta)) + uy*sin(theta)  ]
R = [ uy*ux*(1-cos(theta)) + uz*sin(theta)  cos(theta)+uy^2*(1-cos(theta))          uy*uz*(1-cos(theta)) - ux*sin(theta) ]
    [ uz*ux*(1-cos(theta)) - uy*sin(theta)  uz*uy*(1-cos(theta)) + ux*sin(theta)    cos(theta)+uz^2*(1-cos(theta))        ]
*/


void rotationMatrix(float theta, float* u, float* m) {
	float ux = u[0];
	float uy = u[1];
	float uz = u[2];

	M(0,0) = cos(theta)+ux*ux*(1-cos(theta));
	M(0,1) = ux*uy*(1-cos(theta)) - uz*sin(theta);
	M(0,2) = ux*uz*(1-cos(theta)) + uy*sin(theta) ;

	M(1,0) = uy*ux*(1-cos(theta)) + uz*sin(theta);
	M(1,1) = cos(theta)+uy*uy*(1-cos(theta));
	M(1,2) = uy*uz*(1-cos(theta)) - ux*sin(theta) ;

	M(2,0) = uz*ux*(1-cos(theta)) - uy*sin(theta);
	M(2,1) = uz*uy*(1-cos(theta)) + ux*sin(theta) ;
	M(2,2) = cos(theta)+uz*uz*(1-cos(theta)) ;

}

void matrixMultiply(float* m, float* x, float* r) {
	r[0] = M(0,0)*x[0] + M(0,1)*x[0] + M(0,2)*x[2];
	r[1] = M(1,0)*x[0] + M(1,1)*x[0] + M(1,2)*x[2];
	r[2] = M(2,0)*x[0] + M(2,1)*x[0] + M(2,2)*x[2];
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
	xprime[0] = x[0];
	xprime[1] = x[1];
	xprime[2] = x[2];

}

void transform(float* xprime, float *x) {

    float d_foot = cone.foot;
    float d_rim = cone.mouth;
    float h = cone.height ;
	float w = cone.w;
    float f = d_foot * M_PI / w ; 
    float s = (d_foot*h)/(d_rim-d_foot);
    float r1 = sqrt( s*s + 0.25*d_foot*d_foot );
    float sin_phi = 0.5*d_foot / r1;
    float cos_phi = s / r1;
    //float th = 0.0*f;
    float thicknessf = 1.0;
    float theta0 = d_foot * M_PI / r1;
    float a=r1;
    float c=w/theta0;

	//float max_z = 0;
    //float min_z = 0;

	float xx = x[0];
	float yy = x[1];
	float zz = x[2];

	float z = exp( yy/c)*a-a;
    if (cone.max_z < z - thicknessf*f*zz*sin_phi) {cone.max_z = z - thicknessf*f*zz*sin_phi;}
    if (cone.min_z > z - thicknessf*f*zz*sin_phi) {cone.min_z = z - thicknessf*f*zz*sin_phi;}
    float r = (z/h)*d_rim/2 + (1.0-z/h)*d_foot/2 + thicknessf*f*zz*cos_phi;
    float theta = 2.0*M_PI*( xx / w );

	xprime[0] = (r*cos(theta)); 
	xprime[1] = (r*sin(theta)) ;//
	xprime[2] = z - thicknessf*f*x[2]*sin_phi;
    //print "      vertex", (r*cos(theta)), (r*sin(theta)), z - thicknessf*f*$4*sin_phi

}

//}

float norm(float* n) {
	return sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );
}

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
	n[0] /= l;
	n[1] /= l;
	n[2] /= l;
}


void print_triangle_raw( float* v1, float* v2, float* v3 ) {


	float n[3];
	triangle_normal(v1,v2,v3,&n[0]);

	//printf("  facet normal %f %f %f\n", n[0], n[1], n[2]);
	printf("  facet normal %f %f %f\n", 0.0,0.0,0.0);
	printf("    outer loop\n");
	printf("      vertex %f %f %f\n",v1[0],v1[1],v1[2]);
	printf("      vertex %f %f %f\n",v2[0],v2[1],v2[2]);
	printf("      vertex %f %f %f\n",v3[0],v3[1],v3[2]);
	printf("    endloop\n");
	printf("  endfacet\n");


}



void print_triangle( float* z1, float* z2, float* z3 ) {

	float v1[3];
	float v2[3];
	float v3[3];

#if 1
	transform(v1,z1);
	transform(v2,z2);
	transform(v3,z3);
#else
	idTransform(v1,z1);
	idTransform(v2,z2);
	idTransform(v3,z3);
#endif
	float n[3];
	triangle_normal(v1,v2,v3,&n[0]);

	printf("  facet normal %f %f %f\n", n[0], n[1], n[2]);
	printf("    outer loop\n");
	printf("      vertex %f %f %f\n",v1[0],v1[1],v1[2]);
	printf("      vertex %f %f %f\n",v2[0],v2[1],v2[2]);
	printf("      vertex %f %f %f\n",v3[0],v3[1],v3[2]);
	printf("    endloop\n");
	printf("  endfacet\n");

}

void print_cone_triangles(float thickness) {

    float d_foot = cone.foot;
    float d_rim = cone.mouth;
    float h = cone.height ;
	float w = cone.w;
    float f = d_foot * M_PI / w ; 
    float s = (d_foot*h)/(d_rim-d_foot);
    float r1 = sqrt( s*s + 0.25*d_foot*d_foot );
    //float sin_phi = 0.5*d_foot / r1;
    float cos_phi = s / r1;

    if (cone.max_z < h) {cone.max_z = h;}

    int n_segs = 200; 
    float d_theta = 2.0*M_PI/n_segs;
    for (int i = 0; i < n_segs; ++i) {
        float theta1 = i*d_theta;
        float theta2 = theta1+d_theta;
        float h1 = cone.min_z;
        float h2 = cone.max_z;
        float r1 = (cone.min_z/h)*d_rim/2.0 + (1.0-cone.min_z/h)*d_foot/2.0 - f*thickness/cos_phi ;
        float r2 = (cone.max_z/h)*d_rim/2.0 + (1.0-cone.max_z/h)*d_foot/2.0 - f*thickness/cos_phi ;

		//float p1[3] = { (r1-th)*cos(theta1), (r1-th)*sin(theta1), h1 };
		//float p2[3] = { (r1-th)*cos(theta2),  (r1-th)*sin(theta2), h1 };
		//float p3[3] = { (r2-th)*cos(theta1),  (r2-th)*sin(theta1), h2};

		float c1[3] = {0,0,h1};
		float c2[3] = {0,0,h2};
		float p1[3] = { (r1)*cos(theta1), (r1)*sin(theta1), h1 };
		float p2[3] = { (r1)*cos(theta2),  (r1)*sin(theta2), h1 };
		float p3[3] = { (r2)*cos(theta1),  (r2)*sin(theta1), h2};
		float p4[3] = { (r2)*cos(theta2),  (r2)*sin(theta2), h2};

		print_triangle_raw(&p1[0],&p2[0],&p3[0]);
		print_triangle_raw(&p3[0],&p2[0],&p4[0]);
		print_triangle_raw(&c1[0],&p2[0],&p1[0]);
		print_triangle_raw(&c2[0],&p3[0],&p4[0]);
	}


/*
        
    #outer tall triangles
        print "  facet normal " cos(theta1), sin(theta1), 0
        print "    outer loop"
        print "      vertex ", (r1-th)*cos(theta1),  (r1-th)*sin(theta1), h1
        print "      vertex ", (r1-th)*cos(theta2),  (r1-th)*sin(theta2), h1
        print "      vertex ", (r2-th)*cos(theta1),  (r2-th)*sin(theta1), h2
        print "    endloop"
        print "  endfacet"

        print "  facet normal " cos(theta1), sin(theta1), 0
        print "    outer loop"
        print "      vertex ", (r2-th)*cos(theta2),  (r2-th)*sin(theta2), h2
        print "      vertex ", (r2-th)*cos(theta1),  (r2-th)*sin(theta1), h2
        print "      vertex ", (r1-th)*cos(theta2),  (r1-th)*sin(theta2), h1
        print "    endloop"
        print "  endfacet"

    #top   triangles
        print "  facet normal 0 0 1"
        print "    outer loop"
        print "      vertex ", 0, 0, h2
        print "      vertex ", (r2-th)*cos(theta1), (r2-th)*sin(theta1), h2
        print "      vertex ", (r2-th)*cos(theta2), (r2-th)*sin(theta2), h2
        print "    endloop"
        print "  endfacet"

    #bottom  triangles
        print "  facet normal 0 0 -1"
        print "    outer loop"
        print "      vertex ", 0, 0, h1
        print "      vertex ", (r1-th)*cos(theta2), (r1-th)*sin(theta2), h1
        print "      vertex ", (r1-th)*cos(theta1), (r1-th)*sin(theta1), h1
        print "    endloop"
        print "  endfacet"
*/

}

float dot(float* x, float* y) {
	return ( x[0]*y[0]+x[1]*y[1]+x[2]*y[2] );
}

void print_patch( float* x1,float* x2,float* x3, float thickness, float chamfer) {
	float t1[3] = { x2[0] - x1[0], x2[1] - x1[1] , 0 };
	float t2[3] = { x3[0] - x2[0], x3[1] - x2[1] , 0 };
	float normal[3] = { 0,0, -1 };

	float tn1[3];
	float tn2[3];

	float m1[9];
	float m2[9];

	float t1xt2[3];

	normalize(&t1[0],&tn1[0]);
	normalize(&t2[0],&tn2[0]);

	cross_product(&tn1[0],&tn2[0],&t1xt2[0]);

	if (dot(&t1xt2[0], normal ) <0.0) {return ; }

	rotationMatrix(chamfer, tn1, &m1[0]);
	rotationMatrix(chamfer, tn2, &m2[0]);

	float chamfer_v1[3]; //  = { 0,0, 1 };
	float chamfer_v2[3]; //  = { 0,0, 1 };

	matrixMultiply(&m1[0],&normal[0],&chamfer_v1[0]);
	matrixMultiply(&m2[0],&normal[0],&chamfer_v2[0]);
	float r = thickness / cos(chamfer);

	float v1[3] = { x2[0], x2[1], 0 } ;
	//float v2[3] = { this_x, this_y, 0 } ;
	float v2[3] = { x2[0] + r*chamfer_v1[0], x2[1] + r*chamfer_v1[1], r*chamfer_v1[2]  } ;
	float v3[3] = { x2[0] + r*chamfer_v2[0], x2[1] + r*chamfer_v2[1], r*chamfer_v2[2]  } ;

	print_triangle(v1,v3,v2); 


}


void print_quad(float last_x, float last_y, float this_x, float this_y, float thickness, float chamfer) {

	float t[3] = { this_x - last_x, this_y - last_y , 0 };
	float normal[3] = { 0,0, -1 };
	float chamfer_v[3]; //  = { 0,0, 1 };
	float tn[3] ;
	float m[9];
	normalize(&t[0],&tn[0]);
	rotationMatrix(chamfer, tn, &m[0]);
	matrixMultiply(&m[0],&normal[0],&chamfer_v[0]);

	float r = thickness / cos(chamfer);
	float v1[3] = { last_x, last_y, 0 } ;
	float v2[3] = { this_x, this_y, 0 } ;
	float v3[3] = { last_x + r*chamfer_v[0], last_y + r*chamfer_v[1], r*chamfer_v[2]  } ;
	float v4[3] = { this_x + r*chamfer_v[0], this_y + r*chamfer_v[1], r*chamfer_v[2]  } ;

//	v1 = {}; 

	print_triangle(v1,v2,v3); 
	print_triangle(v3,v2,v4); 


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
	float fb[8]  ;

	float cx,cy;
	//float t = 0.0f, pt = 0.0f;
	TESSalloc ma;
	TESStesselator* tess = 0;
	const int nvp = 3;
	unsigned char* vflags = 0;

	struct MemPool pool;
	unsigned char* mem;   // [1024*1024*20];
	int nvflags = 0;

	//printf("hi %s %s\n",argv[0],argv[1]); fflush(stdout);
	//bg = svgParseFromFile(argv[1]);
	bg = nsvgParseFromFile(argv[1], "px", 1136.0f);

	if (!bg) { printf("error parsing %s\n",argv[1]); return -1; }
	float thickness = atof(argv[2]);
	cone.thickness = thickness;
	chamfer = atof(argv[7]);

	if (argc>3) {
		cone.foot = atof(argv[3]);
		cone.mouth = atof(argv[4]);
		cone.height = atof(argv[5]);
		cone.w = atof(argv[6]);
		cone.max_z =0.0;
		cone.min_z = 0.0;
	}

	//printf("#thickness %f\n",thickness);
	printf("solid x\n");
	//bounds[0] = bounds[2] = bg->pts[0];
	//bounds[1] = bounds[3] = bg->pts[1];
	int np = 0;
	int done_init=0;

	for (NSVGshape *shape = bg->shapes; shape != NULL; shape = shape->next) {
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {
			//printf("x:\nx:\n");
			for (int i = 0; i < path->npts; ++i) {
				float* p = &path->pts[i*2];
				const float x = p[0];
				const float y = p[1];

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
				//drawCubicBez(p[0],p[1], p[2],p[3], p[4],p[5], p[6],p[7]);
			}
		}
	}
/*
	for (it = bg; it != NULL; it = it->next)
		{
			//printf("svg element:\n");

			for (i = 0; i < it->npts; ++i)
			{
				const float x = it->pts[i*2];
				const float y = it->pts[i*2+1];
				//printf("%f %f\n",x,y);
				printf( "x: %f %f\n", x,y );

				np+=1;
				if (x < bounds[0]) bounds[0] = x;
				if (y < bounds[1]) bounds[1] = y;
				if (x > bounds[2]) bounds[2] = x;
				if (y > bounds[3]) bounds[3] = y;
			}
			
		}*/
	cx = (bounds[0]+bounds[2])/2;
	cy = (bounds[1]+bounds[3])/2;
	width = bounds[2]-bounds[0]; 
	height = bounds[3]-bounds[1];

	cone.w = width;
	cone.HH = bounds[3]-bounds[1]-1.0;
	//printf("width:%f\n",cone.w);
	//printf("bounds: %f,%f %f,%f\n",bounds[0],bounds[1],bounds[2],bounds[3]); fflush(stdout);
	//printf("center: %f,%f\n",cx,cy); fflush(stdout);
	//printf("n input vertices: %d\n",np); fflush(stdout);

	mem = malloc( np*4096 );
	//printf("mem=%x\n",mem);
	memset(mem,0, np*4096 );
	//printf("allocate %d\n",np*1024); fflush(stdout);

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
			for (NSVGshape *shape = bg->shapes; shape != NULL; shape = shape->next) {
				for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {
					tessAddContour(tess, 2, path->pts, sizeof(float)*2, path->npts);
				}
			}
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
						//printf("tc:\n");
						//printf("tc:\n");
				//	printf("add contour %d/%d %d %d\n",i,nelems-1,b,n);

/*
					for (j = 0; j < n; ++j)
					{
						printf("tc: %f %f\n",verts[b*2+j*2], verts[b*2+j*2+1]);
					}
					*/
		

					tessAddContour(tess, 2, &verts[b*2], sizeof(float)*2, n);
				}

				
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

				}
			}

			if (tessTesselate(tess, TESS_WINDING_ABS_GEQ_TWO, TESS_BOUNDARY_CONTOURS, 0, 0, 0)) {

				const float* verts = tessGetVertices(tess);
				const int* vinds = tessGetVertexIndices(tess);
				const int nverts = tessGetVertexCount(tess);
				const int* elems = tessGetElements(tess);
				const int nelems = tessGetElementCount(tess);

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
				//	printf("add contour %d/%d %d %d\n",i,nelems-1,b,n);

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

					float last_x = verts[b*2];
					float last_y = verts[b*2+1];

					float this_x; float this_y;
					for (j = 1; j < n; ++j)
					{
						this_x = verts[b*2+j*2];
						this_y = verts[b*2+j*2+1];

						print_quad( last_x, last_y, this_x, this_y, thickness, chamfer);

						last_x = this_x;
						last_y = this_y;
						//printf("tc: %f %f\n",verts[b*2+j*2], verts[b*2+j*2+1]);
					}

					this_x = verts[b*2];
					this_y = verts[b*2+1];
					print_quad( last_x, last_y, this_x, this_y, thickness, chamfer );

					//tessAddContour(tess, 2, &verts[b*2], sizeof(float)*2, n);
				}

				//side wall patches between chamfered rectangles
				for (i = 0; i < nelems; ++i)
				{
				
					int b = elems[i*2];
					int n = elems[i*2+1];

					for (j = 0; j < n; ++j)
					{
						float x1[2] = { verts[b*2+(j%n)*2], verts[b*2+(j%n)*2+1] } ;
						float x2[2] = { verts[b*2+((j+1)%n)*2], verts[b*2+((j+1)%n)*2+1] } ;
						float x3[2] = { verts[b*2+((j+2)%n)*2], verts[b*2+((j+2)%n)*2+1] } ;
						//float x4[2] = { verts[b*2+((j+3)%n)*2], verts[b*2+((j+3)%n)*2+1] } ;

						print_patch( x1,x2,x3, thickness, chamfer);

					}

				}


			//for (it = bg; it != NULL; it = it->next)
			//	tessAddContour(tess, 2, it->pts, sizeof(float)*2, it->npts);
			
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

				}

				//printf("done adding contours\n"); //TESS_WINDING_POSITIVE,
				if (!tessTesselate(tess, TESS_WINDING_ABS_GEQ_TWO, TESS_POLYGONS, nvp, 2, 0)) {
				//if (!tessTesselate(tess,TESS_WINDING_ABS_GEQ_TWO, TESS_POLYGONS, nvp, 2, 0))
					tess = 0;
					printf("tesselate returned zero 2nd time\n");
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

					// Draw polygons.
					//glColor4ub(255,255,255,128);
					// top
					for (i = 0; i < nelems; ++i)
					{

						const int* p = &elems[i*nvp];

						float v1[3] = { verts[p[0]*2],verts[p[0]*2+1],-thickness };
						float v2[3] = { verts[p[1]*2],verts[p[1]*2+1],-thickness };
						float v3[3] = { verts[p[2]*2],verts[p[2]*2+1],-thickness };

						print_triangle(v1,v2,v3);

						v1[2] = 0.0;
						v2[2] = 0.0;
						v3[2] = 0.0;

						print_triangle(v1,v3,v2);

					}

				}

		}
	
	
	print_cone_triangles(thickness);
	
	
	
	printf("endsolid x\n");



	return 0;
}
