#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "meshindex2d.h"
#include "surface.h"
#include "triangles.h"
#include "vectorops.h"


pottoy_spec_t spec;  
void transf_X(float* xyz, float* uv ) {
	faceted_X(xyz,uv, &spec ) ;
};

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

//void clip_triangle( float* a, float* b, float* c, float* d, Plane* p, int* nt ) {

void print_triangle( float* z1, float* z2, float* z3 , FILE* fp) ;

void printf_triangle(float* x) {
    printf("t: %.2f %.2f %.2f ; %.2f %.2f %.2f ; %.2f %.2f %.2f \n", x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8] );
}

#define DEBUG_CLIPPING 0
#define TBUFF_SIZE 3*3*64*2
void print_triangle_clipped( float* z1, float* z2, float* z3 , FILE* fp,
         Plane* rs_planes, int n_rs_planes, 
         Plane* xyz_planes, int n_xyz_planes ,
         float* offset   , float droplevel     
         ) {
    float triangles[ TBUFF_SIZE ];
    int triangle_flags[ TBUFF_SIZE ];
    float triangles2[ TBUFF_SIZE ];
    int nt=1;
    int ti=0;
    int r;
    int altered;
    int nt2=0;
    int ti2=0;
    float zhat[3] = {0,0,1};

    char plane_flags[100];

    for(int pi=0;pi<n_rs_planes;pi++) {
        if (dot( rs_planes[pi].n, &zhat[0] )==1.0) {
            plane_flags[pi]=1;
        }
    }
    

    //float drop_offset[3] = { 0,0,-dropdepth };
    set( &triangles[0], z1 );
    set( &triangles[3], z2 );
    set( &triangles[6], z3 );
    //printf("clip triangle\n "); fflush(stdout);

    memset(&triangle_flags[0],0,sizeof(int)*TBUFF_SIZE);
    while(ti<nt) {
        int nr = 1;
        int nt0 = nt;
        for(int pi=0;pi<n_rs_planes;pi++) {
            int tr = 0; 
        //  printf("clip triangle %d %d %d %f %f %f -- %f %f %f -- %f %f %f\n",ti,nt,triangle_flags[ti],z1[0],z1[1],z1[2],z2[0],z2[1],z2[2],z3[0],z3[1],z3[2]); fflush(stdout);

            if (!triangle_flags[ti]) {
                clip_triangle( &triangles[ ti*9 ] , &triangles[ ti*9 +3], &triangles[ ti*9 +6], 
                &triangles[ nt*9 ], &rs_planes[pi], &tr, &altered );
            } else { tr=1; }
            nr *= tr;


            int n_that_drop=0;


            if (tr>1) {
                // if we're adding new sides perpendicular to the cutting plane...  
                nt++; // //printf("nt=%d\n",nt);
                if ((droplevel!=0.0) && (!plane_flags[pi])) {
                   // int ntt = nt;
                    // for (int i=ti; i<ntt; i++) {
                        //printf("####\n");

                        r=add_side_quad( &triangles[ (nt-1)*9 ] ,  &rs_planes[pi], &triangles[ nt*9 ], droplevel ); if (r>0) {
                            nt += r;
                            for(int j=0; j<r; j++)  {triangle_flags[nt-1-j]=1; }
                            n_that_drop+=2;
                        } else { 
                           // printf("no drops 1\n");
                            }
//                    }
                }
              //  printf("done sidequads %d %d\n",nt,TBUFF_SIZE); fflush(stdout);
                assert(nt*9<TBUFF_SIZE);
             } 
 
            if ((tr>=1)&& altered &&(droplevel!=0.0)&& (!plane_flags[pi])) {
                r=add_side_quad( &triangles[ ti*9     ] ,  &rs_planes[pi], &triangles[ nt*9 ], droplevel ); if (r>0) {
                    nt += r;
                    for(int j=0; j<r; j++)  {triangle_flags[nt-1-j]=1; }
                    n_that_drop+=1;
                }else {
                  //  printf("no drops 0\n");
                    }
            }
      //     printf("xydrops %d\n", n_that_drop);

        }
      //  printf("ti=%d ti*9+6=%d %d; n rs planes: %d\n",ti,ti*9 +6, TBUFF_SIZE, n_rs_planes);fflush(stdout);

        if (nr>0) {
            //print_triangle( &triangles[ ti*9 ] , &triangles[ ti*9 +3], &triangles[ ti*9 +6], fp );
            assert( ti*9 +6 < TBUFF_SIZE );
        //    printf("ti=%d ti*9+6=%d %d\n",ti,ti*9 +6, TBUFF_SIZE);fflush(stdout);
            transform( &triangles2[nt2*9] , &triangles[ ti*9 ]);
            transform( &triangles2[nt2*9+3] , &triangles[ ti*9 +3]);
            transform( &triangles2[nt2*9+6] , &triangles[ ti*9 +6]);

            if (( fabs( triangles2[nt2*9]) > 20000.0 )||
            fabs( triangles2[nt2*9+1]) > 20000.0  ||
            fabs( triangles2[nt2*9+2]) > 20000.0 
            )  {printf( "blowup1: %f %f %f -> %f %f %f  \n"
                , triangles[ti*9], triangles[ti*9+1], triangles[ti*9+2]
                , triangles2[nt2*9], triangles2[nt2*9+1], triangles2[nt2*9+2]
            );  }
            if (( fabs( triangles2[nt2*9+3]) > 20000.0 )||
            fabs( triangles2[nt2*9+4]) > 20000.0  ||
            fabs( triangles2[nt2*9+5]) > 20000.0 
            ) {printf( "blowup1: %f %f %f -> %f %f %f \n"
                , triangles[ti*9+3], triangles[ti*9+4], triangles[ti*9+5]
                , triangles2[nt2*9+3], triangles2[nt2*9+4], triangles2[nt2*9+5]
            );  }
            if (( fabs( triangles2[nt2*9+6]) > 20000.0 )||
            fabs( triangles2[nt2*9+7]) > 20000.0  ||
            fabs( triangles2[nt2*9+8]) > 20000.0 
            ) {printf( "blowup1: %f %f %f -> %f %f %f  \n"
                , triangles[ti*9+6], triangles[ti*9+7], triangles[ti*9+8]
                , triangles2[nt2*9+6], triangles2[nt2*9+7], triangles2[nt2*9+8]
            );  }

            
            
            //}

            nt2++;
        } 

        ti++;
    }
   // printf("n coming out of rs clips %d\n ",ti2); fflush(stdout);

    int n_written=0;
    while(ti2<nt2) {
        int nr = 1;
        int nt20 = nt2;
        for(int pi=0;pi<n_xyz_planes;pi++) {
            //printf("clip with xyzplane ");
           // printf_triangle(&triangles2[ ti2*9 ]);
            int tr = 0; 
            clip_triangle( &triangles2[ ti2*9 ] , &triangles2[ ti2*9 +3], &triangles2[ ti2*9 +6], 
            &triangles2[ nt2*9 ], &xyz_planes[pi], &tr, &altered );
            nr *= tr;
            if (tr>1) { 
              //  printf_triangle(&triangles2[ nt2*9 ]);
            nt2++; //printf("nt=%d\n",nt);
             assert(nt*9<TBUFF_SIZE);} 
        }
        if (nr>0) {
                add( &triangles2[ ti2*9 ], &triangles2[ ti2*9 ],  offset);
                add( &triangles2[ ti2*9+3 ], &triangles2[ ti2*9+3 ],  offset);
                add( &triangles2[ ti2*9+6 ], &triangles2[ ti2*9+6 ],  offset);

            stl_triangle( 
                &triangles2[ ti2*9 ]  , 
                &triangles2[ ti2*9 +3] , 
                &triangles2[ ti2*9 +6] , fp );
                fflush(fp);
            n_written++;
        } 

        ti2++;
    }

 //   printf("done clip triangle %d\n ",n_written); fflush(stdout);


    //printf("       done clip triangle\n ");

}

void print_triangle_bbox( float* z1, float* z2, float* z3 , FILE* fp, float* bbox) {

	float v1[3];
	float v2[3];
	float v3[3];

//    printf("transforming %f %f %f\n", z1[0],z1[1],z1[2]);
//    printf("transforming %f %f %f\n", z2[0],z2[1],z2[2]);
//    printf("transforming %f %f %f\n", z3[0],z3[1],z3[2]);
//    fflush(stdout);
	transform(v1,z1);
	transform(v2,z2);
	transform(v3,z3);
 //   printf("     %f %f %f\n", v1[0],v1[1],v1[2]);
 //   printf("     %f %f %f\n", v2[0],v2[1],v2[2]);
 //   printf("     %f %f %f\n", v3[0],v3[1],v3[2]);


   // fflush(stdout);

    float miny = bbox[1];
    float maxy = bbox[4];

    if (v1[1] > maxy && v2[1] > maxy && v3[1] > maxy  ) {return;}
    if (v1[1] < miny && v2[1] < miny  && v3[1] < miny ) {return;}
    if (isnan( v1[0] )) {return;}
    if (isnan( v1[1] )) {return;}
    if (isnan( v1[2] )) {return;}
    if (isnan( v2[0] )) {return;}
    if (isnan( v2[1] )) {return;}
    if (isnan( v2[2] )) {return;}
    if (isnan( v3[0] )) {return;}
    if (isnan( v3[1] )) {return;}
    if (isnan( v3[2] )) {return;}

	float n[3];
	triangle_normal(v1,v2,v3,&n[0]);
//    printf("print triangle with normal %f %f %f\n",n[0], n[1], n[2] );
	if (n[0]==0.0 && n[1]==0.0 && n[2]==0.0) { 
        //printf("zero normla\n") ;; 
        return;}

	fprintf(fp,"  facet normal %f %f %f\n", n[0], n[1], n[2]);
	fprintf(fp,"    outer loop\n");
	fprintf(fp,"      vertex %f %f %f\n",v1[0],-v1[2],v1[1]);
	fprintf(fp,"      vertex %f %f %f\n",v2[0],-v2[2],v2[1]);
	fprintf(fp,"      vertex %f %f %f\n",v3[0],-v3[2],v3[1]);
	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");
    fflush(fp);

}
/*

void print_triangle( float* z1, float* z2, float* z3 , FILE* fp) {

	float v1[3];
	float v2[3];
	float v3[3];

	transform(v1,z1);
	transform(v2,z2);
	transform(v3,z3);

	float n[3];
	triangle_normal(v1,v2,v3,&n[0]);
//    printf("print triangle with normal %f %f %f\n",n[0], n[1], n[2] );
	if (n[0]==0.0 && n[1]==0.0 && n[2]==0.0) { 
        //printf("zero normla\n") ;; 
        return;}

	fprintf(fp,"  facet normal %f %f %f\n", n[0], n[1], n[2]);
	fprintf(fp,"    outer loop\n");
	fprintf(fp,"      vertex %f %f %f\n",v1[0],-v1[2],v1[1]);
	fprintf(fp,"      vertex %f %f %f\n",v2[0],-v2[2],v2[1]);
	fprintf(fp,"      vertex %f %f %f\n",v3[0],-v3[2],v3[1]);

//	fprintf(fp,"      vertex %f %f %f\n",v1[0],v1[1],v1[2]);
//	fprintf(fp,"      vertex %f %f %f\n",v2[0],v2[1],v2[2]);
//	fprintf(fp,"      vertex %f %f %f\n",v3[0],v3[1],v3[2]);
	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");

}*/


int main(int argc, char *argv[])
{

	int n_sectors = 90;
	int n_levels = 50;
    int mold_mode = 0;


    if (!( argc==4 )) {
        //printf("usage: geowrap <input.json> <output.obj>\n");
        printf("usage: meshwrap <texture_geometry.obj> <surface.obj> <out.stl>\n"); //<input.json> <output.obj>\n");
        exit(0);
    }

    if (argc==4) {

        // argv[1] -- .obj file to be wrapped, or .svg file to be embossed
        // argv[2] -- obj file of mesh to wrap around with the uv coordinates from lscm
        // argv[3] -- output stl file
 
        FILE* stlfile = fopen(argv[3], "wt"); 
//	    fprintf(stlfile,"solid x\n");

        MeshTriangles* texture = parse_triangles_with_normals(argv[1],(float) 1.0);
        if (texture == NULL) {return 1;} 

        float width = mesh_width(texture);

        MeshTriangles* mt = parse_triangles_with_normals(argv[2],width);
        if (mt == NULL) {return 1;} 
        meshindex* mi = build_mesh_index( mt->points, mt->triangles, mt->ntriangles, 10, 10 ) ;

        printf("parsed\n");
//        return(0);

        for (int i=0;i<texture->nxpoints;i++) {
            int hc=0;
           mesh_interpolation_xyz(mt,&texture->xpoints[i*3], &texture->xpoints[i*3],mi,&hc );    
           if (0==i%10000) {printf("p %d/%d  %f\n",i,texture->nxpoints, (((float ) i)/((float) texture->nxpoints)));}        
        }


        write_to_obj( texture, stlfile, 0) ;


        fclose(stlfile);
        free_meshtriangles( mt ) ;

        free_meshtriangles( texture ) ;
    }
    printf("done ok\n");
    return(0);
}

