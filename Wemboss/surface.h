#include "tesselator.h"



void X(float* xyz, float* vtheta);
//void Xs(float* xyz, float* uv, pottoy_spec_t* spec ) ;


int write_surface_obj(char* fn);
//int write_surface_obj_spec(spec,argv[2]);   


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

static void scan_array(const char *str, int len, void *user_data);
int read_spec(char* filename, pottoy_spec_t* spec );

void faceted_X(float* xyz, float* uv, pottoy_spec_t* spec ) ;

void surface_obj_bbox(float* bbox, pottoy_spec_t* spec) ;
void scale_spec(float scale, pottoy_spec_t* spec) ;


void add_triangle_contours(TESStesselator* tess, pottoy_spec_t* spec, int mn_sectors, int mn_levels);

void mold_parting_norm_ray( float u, float l,float thickness, int quadrant  , pottoy_spec_t* spec, float* result) ;

//void mold_parting_norm_ray( float u, float l, int quadrant  , pottoy_spec_t* spec, float* result) ;
void write_parting_sufrace_stl( int quadrant , float l, pottoy_spec_t* spec, FILE* stlfile, float* offset, 
     float thickness,  int reverse) ;

//void write_parting_sufrace_stl( int quadrant, float l , pottoy_spec_t* spec, FILE* stlfile) ;
void stl_triangle( float* v1, float* v2, float* v3, FILE* fp ) ;

void write_surface_stl( void(*trnsfrm)(float*, float*), FILE* stlfile,float mins, float maxs, float thickness , float* offset) ;
