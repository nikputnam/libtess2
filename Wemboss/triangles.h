#include <stdbool.h>

#ifndef _tringles_h_
#define _tringles_h_


#include "meshindex2d.h"

typedef struct MeshTriangles
{
	int ntriangles;				
	int npaths;				
	int npoints;				
	int nxpoints;				
	int ntpoints;				
	float* points;			
	float* xpoints;			
	float* tpoints;			
	int* triangles;		
	int* ttriangles;		
    float* paths;	

	int max_ntriangles;				
	int max_npaths;				
	int max_npoints;				
	int max_nxpoints;				
	int max_ntpoints;				

} MeshTriangles;



typedef struct Plane {
	float n[3];
	float d;
} Plane ; 


MeshTriangles* parse_triangles(char* filename, float width) ;
MeshTriangles* parse_triangles_internal(char* filename, float width,int with_normals, int reduplicate) ;
MeshTriangles* parse_triangles_with_normals(char* filename, float width) ;
void mesh_interpolation(MeshTriangles* mt, float* xy, float* uv, meshindex* mi, int* hit_count) ;
void write_to_stl( MeshTriangles* t, FILE* stlfile ) ;
void write_to_obj(MeshTriangles* t, FILE* stlfile) ;

void add_triangle(float* v1, float* v2, float* v3 , MeshTriangles* mt);

void print_triangle_raw( float* v1, float* v2, float* v3 , FILE* fp) ;
void triangle_normal(float* v1, float* v2, float* v3, float *n) ;

void apply_interpolation_to_mesh( MeshTriangles* texture, MeshTriangles* mt , float zscale,meshindex* mi);
void apply_transform_to_mesh(MeshTriangles*  texture,  void(*trnsfrm)(float*, float*) );

float signed_distance_to_plane(float* a, Plane* p);
bool segment_spans_plane( float* a, float* b, Plane* p );

void clip_triangle( float* a, float* b, float* c, float* d, Plane* p, int* nt, int* altered ) ;

bool triangle_spans_plane( float* a, float* b, float* c, Plane* p ) ;


void segment_plane_intersection( float* a, float* b, float* q, Plane* p ) ;

int add_side_quad( float* triangle ,  Plane* p, float* new_triangles, float droplevel );
void free_meshtriangles( MeshTriangles* mt ) ;

#endif