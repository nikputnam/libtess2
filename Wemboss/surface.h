#include "tesselator.h"
#include "triangles.h"

#ifndef _surface_h_
#define _surface_h_


void X(float* xyz, float* vtheta);
//void Xs(float* xyz, float* uv, pottoy_spec_t* spec ) ;


int write_surface_obj(char* fn);
//int write_surface_obj_spec(spec,argv[2]);   

typedef struct stl_output_config {
    FILE* fp;
    float scaling;
    Plane* clipping_planes;
    int n_clipping_planes;
} stl_output_config;


typedef struct pottoy_struct {

    int facet;
    int n_facets;
    float puff;
    float squash;
    float twist;
    int npoints;
    float* points;
/*
    "facet":true,
    "n_facets":9,
    "points":[{"x":289.76843686945404,"y":0},{"x":414.91846107193595,"y":191.8845958690277},
        {"x":626.0426400209401,"y":388.4019456777616},{"x":445.835534762008,"y":698.9204108718246}],
    "puff":-1.8,
    "squash":1,
    "twist":0.17,
    "npoints":4
*/

} pottoy_spec_t ;

void Xs(float* xyz, float* vtheta, pottoy_spec_t* spec );
//int write_surface_obj2(char* fn, pottoy_spec_t* spec) ;
int write_surface_obj2(char* fn, pottoy_spec_t* spec, int n_sectors, int n_levels) ;
int write_surface_obj0(char* fn, pottoy_spec_t* spec, int n_sectors, int n_levels) ;

static void scan_array(const char *str, int len, void *user_data);
int read_spec(char* filename, pottoy_spec_t* spec );

void faceted_X(float* xyz, float* uv, pottoy_spec_t* spec ) ;

void surface_obj_bbox(float* bbox, pottoy_spec_t* spec) ;
void scale_spec(float scale, pottoy_spec_t* spec) ;


void add_triangle_contours(TESStesselator* tess, pottoy_spec_t* spec, int mn_sectors, int mn_levels);

void mold_parting_norm_ray( float u, float l,float thickness, int quadrant  , pottoy_spec_t* spec, float* result) ;

//void mold_parting_norm_ray( float u, float l, int quadrant  , pottoy_spec_t* spec, float* result) ;
void write_parting_sufrace_stl( int quadrant , float l, pottoy_spec_t* spec, stl_output_config cf, float* offset, 
     float thickness,  int reverse) ;

//void write_parting_sufrace_stl( int quadrant, float l , pottoy_spec_t* spec, FILE* stlfile) ;
void stl_triangle( float* v1, float* v2, float* v3, FILE* fp ) ;
void output_triangle( float* v1, float* v2, float* v3, stl_output_config cf ) ;

void write_texture_back_stl( void(*trnsfrm)(float*, float*), 
            float mins, float maxs , float thickness, float* offset,
            stl_output_config cf
            ) ;

void write_surface_stl( void(*trnsfrm)(float*, float*), FILE* stlfile,float mins, float maxs, float thickness , float* offset) ;
void write_floor_flange_stl(  pottoy_spec_t* spec, FILE* stlfile,float mins, float maxs ,float length, float thickness, float* offset,
 float t, int top) ;

int contour_to_mesa(int n, float* p, float* offset, float f, int nt, float* result_buffer) ;


void write_support_ties_stl(
    pottoy_spec_t* spec,
    int quadrant1, int quadrant2, float* offsets1, 
        float* offsets2, float length, float thickness, stl_output_config cf ) ;//
       // &transform, stlfile, mins, maxs, thickness,&offsets[quadrant*3]clup

void output_split_quad( float* ray1 , float* ray2 ,stl_output_config cf , int flip, float a, float b, float c, float d ) ;


#endif