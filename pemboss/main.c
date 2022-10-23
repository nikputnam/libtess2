#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "nanosvg.h"
#include "bez.h"
#include "tess.h"
#include "tesselator.h"
#include "vectorops.h"
#include "triangles.h"

//#include "meshindex2d.h"
//#include "surface.h"
//#include "triangles.h"
//#include "vectorops.h"

void add_triangle_contours0(TESStesselator* tess, int mn_sectors, int mn_levels, float width, float height) {

//return;
    float x[6];
    float uv[2];

#define VV(ii) (height*((double)(ii))/((double) mn_levels))
#define UU(jj) (width*((double) (jj))/((double) (mn_sectors-1)))

        
    for (int j =0; j<mn_levels;j++) {
        for (int i =0; i<mn_sectors;i++) {

            x[0] = UU(i);   x[1] = VV(j);
            x[4] = UU(i+1); x[5] = VV(j);
            x[2] = UU(i);   x[3] = VV(j+1);

            //printf("z %f,%f %f,%f %f,%f \n",x[0],x[1],x[2],x[3],x[4],x[5]);
            tessAddContour(tess, 2, (void*) &(x[0]), sizeof(float)*2, 3);

            x[0] = UU(i)  ; x[1] = VV(j+1);
            x[4] = UU(i+1); x[5] = VV(j);
            x[2] = UU(i+1); x[3] = VV(j+1);

            //printf("z %f,%f %f,%f %f,%f \n",x[0],x[1],x[2],x[3],x[4],x[5]);
            tessAddContour(tess, 2, (void*) &(x[0]), sizeof(float)*2, 3);
        }
    }

 //   fflush(stdout);
}


void add_surface( int mn_sectors, int mn_levels, float width, float height, MeshTriangles* mt, float z, int flip) {

//return;
    float x[9];
//    float uv[2];

#define VV(ii) (height*((double)(ii))/((double) mn_levels))
#define UU(jj) (width*((double) (jj))/((double) (mn_sectors-1)))

        
    for (int j =0; j<mn_levels;j++) {
        for (int i =0; i<mn_sectors;i++) {


            x[0] = UU(i);   x[1] = VV(j);  x[2] = z;
            x[6] = UU(i+1); x[7] = VV(j);   x[8]=z;
            x[3] = UU(i);   x[4] = VV(j+1); x[5] = z;

            //printf("z %f,%f %f,%f %f,%f \n",x[0],x[1],x[2],x[3],x[4],x[5]);
            //tessAddContour(tess, 2, (void*) &(x[0]), sizeof(float)*2, 3);

            if (flip ) {
                add_triangle(&x[0], &x[6], &x[3], mt );
                
            } else {
                add_triangle(&x[0], &x[3], &x[6], mt );
            }

            x[0] = UU(i)  ; x[1] = VV(j+1); x[2] = z;
            x[6] = UU(i+1); x[7] = VV(j);x[8]=z;
            x[3] = UU(i+1); x[4] = VV(j+1);x[5] = z;


            if (flip ) {
            add_triangle(&x[0], &x[6], &x[3], mt );
            } else {
            add_triangle(&x[0], &x[3], &x[6], mt );

            }

            //printf("z %f,%f %f,%f %f,%f \n",x[0],x[1],x[2],x[3],x[4],x[5]);
//            tessAddContour(tess, 2, (void*) &(x[0]), sizeof(float)*2, 3);
        }
    }

 //   fflush(stdout);
}



void add_wedge(float* x20, float* x10,float thickness, float chamfer, MeshTriangles* mt) {
//	return;
    //printf("x10: %f %f %f ; x20: %f %f %f \n",x10[0],x10[1],x10[2],x20[0],x20[1],x20[2]);
    
	float x1[3];// = { x10[0], x10[1], 0 };
	float x2[3];// = { x20[0], x20[1], 0 };
    float xa[3];


	set(&x1[0],x10);
	set(&x2[0],x20);

    x1[2] += thickness ;
    x2[2] += thickness ;

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

	float normal1[3] = { 0,0, 1 };
	float normal2[3] = { 0,0, 1 };
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

	float r = 0.95*thickness / cos(chamfer);
    //r = thickness;
	//float v1[3] = { last_x, last_y, 0 } ;
	//float v2[3] = { this_x, this_y, 0 } ;

	float x3[3] = { x1[0] - r*chamfer_v1[0], x1[1] - r*chamfer_v1[1], x1[2] - r*chamfer_v1[2]  } ;
	float x4[3] = { x2[0] - r*chamfer_v2[0], x2[1] - r*chamfer_v2[1], x2[2] - r*chamfer_v2[2]  } ;

	float x5[3] = { x1[0] - thickness*normal1[0] , x1[1] - thickness*normal1[1] ,x1[2]  - thickness*normal1[2] };
	float x6[3] = { x2[0] - thickness*normal2[0] , x2[1] - thickness*normal2[1] ,x2[2]  - thickness*normal2[2] };

//	v1 = {}; 

	add_triangle(x1,x2,x3, mt); 
	add_triangle(x3,x2,x4, mt); 

	//print_triangle(x1,x3,x5); 
	//print_triangle(x4,x2,x6); 
    
	add_triangle(x6,x5,x3, mt); 
	add_triangle(x4,x6,x3, mt); 
	add_triangle(x2,x1,x5, mt); 
	add_triangle(x2,x5,x6, mt); 
    

}

void add_quad(float last_x, float last_y, float this_x, float this_y, float thickness, MeshTriangles* mt) {

//    return;
	float v1[3] = { last_x, last_y, 0.05*thickness } ;
	float v2[3] = { this_x, this_y, 0.05*thickness } ;
	float v3[3] = { last_x , last_y , thickness  } ;
	float v4[3] = { this_x , this_y , thickness  } ;


	add_triangle(v1,v3,v2, mt); 
	add_triangle(v3,v4,v2, mt); 

}
void add_patch( float* x30,float* x20,float* x10, float thickness, float chamfer, MeshTriangles* mt) {

//    return;
	float x1[3];// = { x10[0], x10[1], 0 };
	float x2[3];// = { x20[0], x20[1], 0 };
	float x3[3];// = { x20[0], x20[1], 0 };
    float xa[3];
    float xb[3];

	if (vequal2(x10,x20)) {printf("x10==x20 nan\n");}
	if (vequal2(x20,x30)) {printf("x20==x30 nan\n");}

	set(&x1[0],x10);
	set(&x2[0],x20);
	set(&x3[0],x30);

    x1[2] += thickness ;
    x2[2] += thickness ;
    x3[2] += thickness ;

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

	float normal2[3] = { 0,0, 1 };


//	float normal2[3]; // = { 0,0, -1 };
//	float normalB[3]; // = { 0,0, -1 };

 //   surface_norm(&normal2[0], &x20[0], transf_X );
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

	float r = 0.95*thickness / cos(chamfer);
    //printf("r: %f \t %f \t %f\n",r,thickness,chamfer);
	float v1[3] = { x2[0], x2[1], x2[2] } ;
	//float v2[3] = { this_x, this_y, 0 } ;
	float v2[3] = { x2[0] - r*chamfer_v1[0], x2[1] - r*chamfer_v1[1], x2[2] - r*chamfer_v1[2]  } ;
	float v3[3] = { x2[0] - r*chamfer_v2[0], x2[1] - r*chamfer_v2[1], x2[2] - r*chamfer_v2[2]  } ;
	float v4[3] = { x2[0] - thickness*normal2[0], x2[1] - thickness*normal2[1], x2[2] - thickness*normal2[2] } ;

	add_triangle(v1,v3,v2, mt); 
	add_triangle(v2,v3,v4, mt); 


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

int main(int argc, char *argv[]) {
    FILE* stlfile;
    float thickness = 5.0;
    NSVGimage *bg = NULL;
	bg = nsvgParseFromFile(argv[1], "px", 98.0f);
	if (!bg) { printf("error parsing %s\n",argv[1]); return -1; }

    	int np = 0;
	int done_init=0;

	float *p0, *p1, *p2, *p3; 

	int n_points = 0;
	int n_paths = 0;
    int ppi = 0;

    float bounds[4];
	float chamfer  = ( (M_PI* 60.0 )/180.0); 
	//coneParams cone;

	float width,height;
    int i,j;
	//float fb[8]  ;
    	float cx,cy;


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


    for (int i=0;i<n_points;i++ ) {
        path_points[2*i+1] = height - path_points[2*i+1];
    }

    printf("bounts %f %f %f %f \n",bounds[0],bounds[1],bounds[2],bounds[3]);
    printf("width: %f ; height: %f\n",width,height);


    // now triangulate
    if (1) {

    MeshTriangles mt;
    mt.max_ntriangles = 1000;
    mt.max_nxpoints = 1000;
    mt.xpoints = (float*) malloc( sizeof(float)*3*mt.max_nxpoints );
    mt.triangles = (int*)   malloc( sizeof(int)*3*mt.max_ntriangles );
    mt.ntriangles =0;
    mt.nxpoints = 0;

    FILE* stlfile = fopen(argv[2], "wt"); 
	fprintf(stlfile,"solid x\n");

    int n_sectors = 90;
    int n_levels = 50;

    TESSalloc ma;
    TESStesselator* tess = 0;
    const int nvp = 3;
    unsigned char* vflags = 0;

    struct MemPool pool;
    unsigned char* mem;   // [1024*1024*20];
    int nvflags = 0;

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
	
#define REVERSE_PATTERN_CONTOURS 1
#if REVERSE_PATTERN_CONTOURS
                tessSetOption(tess,TESS_REVERSE_CONTOURS,1);
#endif

    	for (int i =0; i<n_paths; i++) {
			tessAddContour(tess, 2, &path_points[ path_offsets[i] ] , sizeof(float)*2, path_lengths[i]);
		}
#if REVERSE_PATTERN_CONTOURS
                tessSetOption(tess,TESS_REVERSE_CONTOURS,0);
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

				add_triangle_contours0(tess, n_sectors, n_levels, width, height) ;
				
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
							add_wedge(&x2[0],&x1[0],1.0*thickness,chamfer,&mt);
						} 
#endif					
						add_quad( x1[0], x1[1], x2[0], x2[1], 1.0*thickness, &mt);
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
							add_patch( x3,x2,x1,thickness, chamfer, &mt);
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

				add_triangle_contours0(tess,  n_sectors, n_levels, width, height) ;
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

						float v1[3] = { verts[p[0]*2],verts[p[0]*2+1],1.0*thickness };
						float v2[3] = { verts[p[1]*2],verts[p[1]*2+1],1.0*thickness };
						float v3[3] = { verts[p[2]*2],verts[p[2]*2+1],1.0*thickness };

//						print_triangle(v1,v2,v3, stlfile);
                        add_triangle(v1,v3,v2,&mt);

						v1[2] = 0.05*thickness;
						v2[2] = 0.05*thickness;
						v3[2] = 0.05*thickness;

//						print_triangle(v1,v3,v2, stlfile);
                        add_triangle(v1,v2,v3,&mt);

					}
#endif



/*

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
					*/
				}
		}
	
      add_surface( n_sectors, n_levels, width, height,&mt,  0.1*thickness, 1) ;
      add_surface( n_sectors, n_levels, width, height,&mt, -0.3*thickness, 0) ;

#define YUP 1
    write_to_obj(&mt,stlfile,YUP);


    }




}