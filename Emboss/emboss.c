
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

void transform(float* xprime, float *x) {

    float d_foot = cone.foot;
    float d_rim = cone.mouth;
    float h = cone.height ;
	float w = cone.w;
    float f = d_foot * 3.1415 / w ; 
    float s = (d_foot*h)/(d_rim-d_foot);
    float r1 = sqrt( s*s + 0.25*d_foot*d_foot );
    float sin_phi = 0.5*d_foot / r1;
    float cos_phi = s / r1;
    //float th = 0.0*f;
    float thicknessf = 1.0;
    float theta0 = d_foot * 3.1415 / r1;
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
    float theta = 2.0*3.1415*( xx / w );

	xprime[0] = (r*cos(theta)); 
	xprime[1] = (r*sin(theta)) ;//
	xprime[2] = z - thicknessf*f*x[2]*sin_phi;
    //print "      vertex", (r*cos(theta)), (r*sin(theta)), z - thicknessf*f*$4*sin_phi

}

//}

void triangle_normal(float* v1, float* v2, float* v3, float *n) {

	float a[3];
	float b[3];

	a[0] = v2[0] - v1[0];
	a[1] = v2[1] - v1[1];
	a[2] = v2[2] - v1[2];

	b[0] = v3[0] - v1[0];
	b[1] = v3[1] - v1[1];
	b[2] = v3[2] - v1[2];

	n[0] = a[1]*b[2] - a[2]*b[1];
	n[1] = a[2]*b[0] - a[0]*b[2];
	n[2] = a[0]*b[1] - a[1]*b[0];
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

	transform(v1,z1);
	transform(v2,z2);
	transform(v3,z3);

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

void print_cone_triangles() {
	float h = cone.height;
	//float th = 0.1 ; //; cone.thickness;
    float d_foot = cone.foot;
    float d_rim = cone.mouth;

    if (cone.max_z < h) {cone.max_z = h;}

    int n_segs = 200; 
    float d_theta = 2.0*3.1415/n_segs;
    for (int i = 0; i < n_segs; ++i) {
        float theta1 = i*d_theta;
        float theta2 = theta1+d_theta;
        float h1 = cone.min_z;
        float h2 = cone.max_z;
        float r1 = (cone.min_z/h)*d_rim/2.0 + (1.0-cone.min_z/h)*d_foot/2.0 ;
        float r2 = (cone.max_z/h)*d_rim/2.0 + (1.0-cone.max_z/h)*d_foot/2.0 ;

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
void print_quad(float last_x, float last_y, float this_x, float this_y, float thickness) {

	float v1[3] = { last_x, last_y, thickness } ;
	float v2[3] = { this_x, this_y, thickness } ;
	float v3[3] = { last_x, last_y, -thickness } ;
	float v4[3] = { this_x, this_y, -thickness } ;

//	v1 = {}; 

	print_triangle(v1,v2,v3); 
	print_triangle(v3,v2,v4); 


}

int main(int argc, char *argv[])
{
	struct SVGPath* bg;
	struct SVGPath* it;
	float bounds[4];

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
	bg = svgParseFromFile(argv[1]);
	if (!bg) { printf("error parsing %s\n",argv[1]); return -1; }
	float thickness = atof(argv[2]);
	cone.thickness = thickness;

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
	bounds[0] = bounds[2] = bg->pts[0];
	bounds[1] = bounds[3] = bg->pts[1];
	int np = 0;
	for (it = bg; it != NULL; it = it->next)
		{
			//printf("svg element:\n");

			for (i = 0; i < it->npts; ++i)
			{
				const float x = it->pts[i*2];
				const float y = it->pts[i*2+1];
				//printf("%f %f\n",x,y);
				np+=1;
				if (x < bounds[0]) bounds[0] = x;
				if (y < bounds[1]) bounds[1] = y;
				if (x > bounds[2]) bounds[2] = x;
				if (y > bounds[3]) bounds[3] = y;
			}
			
		}
	cx = (bounds[0]+bounds[2])/2;
	cy = (bounds[1]+bounds[3])/2;
	width = bounds[2]-bounds[0]; 
	height = bounds[3]-bounds[1];

	//cone.w = bounds[2]-bounds[0];
	cone.HH = bounds[3]-bounds[1];
	printf("width:%f\n",cone.w);
	//printf("center: %f,%f\n",cx,cy); fflush(stdout);
	//printf("n input vertices: %d\n",np); fflush(stdout);

	mem = malloc( np*2048 );
	//printf("mem=%x\n",mem);
	memset(mem,0, np*2048 );
	//printf("allocate %d\n",np*1024); fflush(stdout);

	pool.size = 0;
	pool.cap = np*2048; //sizeof(mem);
	pool.buf = mem;
	memset(&ma, 0, sizeof(ma));
	//memset(&ma,0, sizeof(ma))
	ma.memalloc = poolAlloc;
	ma.memfree = poolFree;
	ma.userData = (void*)&pool;
	ma.extraVertices = 2048; // realloc not provided, allow 256 extra vertices.


		pool.size = 0; // reset pool
		tess = tessNewTess(&ma);
		if (tess)
		{
			tessSetOption(tess, TESS_CONSTRAINED_DELAUNAY_TRIANGULATION, 0);

			for (it = bg; it != NULL; it = it->next)
				tessAddContour(tess, 2, it->pts, sizeof(float)*2, it->npts);
			
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

#if 1
		// First combine contours and then triangulate, this removes unnecessary inner vertices.
			if (tessTesselate(tess, TESS_WINDING_ABS_GEQ_TWO, TESS_BOUNDARY_CONTOURS, 0, 0, 0))
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

					//for (j = 0; j < n; ++j)
					//{
					//	printf("tc: %f %f\n",verts[b*2+j*2], verts[b*2+j*2+1]);
					//}

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

						print_quad( last_x, last_y, this_x, this_y, thickness );
						//printf("  facet normal 0 0 1\n");
    					//printf("    outer loop\n");
						//	printf("      vertex %f %f %f\n",verts[p[j]*2],verts[p[j]*2+1],1);
						//printf("    endloop\n");
						//printf("  endfacet\n");

						last_x = this_x;
						last_y = this_y;
						//printf("tc: %f %f\n",verts[b*2+j*2], verts[b*2+j*2+1]);
					}

					this_x = verts[b*2];
					this_y = verts[b*2+1];
					print_quad( last_x, last_y, this_x, this_y, thickness );

					//tessAddContour(tess, 2, &verts[b*2], sizeof(float)*2, n);
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

						v1[2]*=-1.0;
						v2[2]*=-1.0;
						v3[2]*=-1.0;

						print_triangle(v1,v3,v2);

					}

				}

		}
	
	
	print_cone_triangles();
	
	
	
	printf("endsolid x\n");



	return 0;
}
