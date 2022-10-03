#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "surface.h"
#include "triangles.h"
#include "vectorops.h"
#include "meshindex2d.h"


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
          printf("clip triangle %d %d %d %f %f %f -- %f %f %f -- %f %f %f\n",ti,nt,triangle_flags[ti],z1[0],z1[1],z1[2],z2[0],z2[1],z2[2],z3[0],z3[1],z3[2]); fflush(stdout);

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
           printf("xydrops %d\n", n_that_drop);

        }
        printf("ti=%d ti*9+6=%d %d; n rs planes: %d\n",ti,ti*9 +6, TBUFF_SIZE, n_rs_planes);fflush(stdout);

        if (nr>0) {
            //print_triangle( &triangles[ ti*9 ] , &triangles[ ti*9 +3], &triangles[ ti*9 +6], fp );
            assert( ti*9 +6 < TBUFF_SIZE );
            printf("ti=%d ti*9+6=%d %d\n",ti,ti*9 +6, TBUFF_SIZE);fflush(stdout);
            transform( &triangles2[nt2*9] , &triangles[ ti*9 ]);
            transform( &triangles2[nt2*9+3] , &triangles[ ti*9 +3]);
            transform( &triangles2[nt2*9+6] , &triangles[ ti*9 +6]);
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
	fprintf(fp,"      vertex %f %f %f\n",0.1*v1[0],-0.1*v1[2],0.1*v1[1]);
	fprintf(fp,"      vertex %f %f %f\n",0.1*v2[0],-0.1*v2[2],0.1*v2[1]);
	fprintf(fp,"      vertex %f %f %f\n",0.1*v3[0],-0.1*v3[2],0.1*v3[1]);
	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");
    fflush(fp);

}


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
	fprintf(fp,"      vertex %f %f %f\n",0.1*v1[0],-0.1*v1[2],0.1*v1[1]);
	fprintf(fp,"      vertex %f %f %f\n",0.1*v2[0],-0.1*v2[2],0.1*v2[1]);
	fprintf(fp,"      vertex %f %f %f\n",0.1*v3[0],-0.1*v3[2],0.1*v3[1]);

//	fprintf(fp,"      vertex %f %f %f\n",v1[0],v1[1],v1[2]);
//	fprintf(fp,"      vertex %f %f %f\n",v2[0],v2[1],v2[2]);
//	fprintf(fp,"      vertex %f %f %f\n",v3[0],v3[1],v3[2]);
	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");

}


int main(int argc, char *argv[])
{

	int n_sectors = 90;
	int n_levels = 50;
    int mold_mode = 0;


    if (!( argc==3 || argc==6 || argc==7)) {
        printf("usage: geowrap <input.json> <output.obj>\n");
        printf("usage: geowrap <texture_geometry.obj> <input-lscm.obj> <out.stl> <surf.json> <width> [1 to print mold]\n"); //<input.json> <output.obj>\n");
        exit(0);
    }

    if (argc==3) {
		printf("read spec & write obj.\n");
		read_spec(argv[1], &spec );

        if (0) {

            float a[3] = {1,1,0};
            float b[3] = {4,2,0};
            float c[3] = {2,4,0};

            Plane clip_planes[2];
            clip_planes[0] = (Plane){ {0,-1,0}, 3.0 };
            clip_planes[1] = (Plane){ {-1,0,0}, 1.0 };

            //print_triangle_clipped( a,b,c, stdout, clip_planes, 1);

            exit(0);

        }

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

    if (argc>=6) {

        // argv[1] -- .obj file to be wrapped
        // argv[2] -- obj file with the uv coordinates from lscm
        // argv[3] -- output stl file
        // argv[4] -- surface spec in .json
        // argv[5] -- wrapping width of the texture
        // argv[6] (optional) 4-sided mold version
 
        FILE* stlfile = fopen(argv[3], "wt"); 
	    fprintf(stlfile,"solid x\n");

		read_spec(argv[4], &spec );
		float bbox[6];
		surface_obj_bbox(&bbox[0], &spec) ;
		printf("bbox: x:  %f - %f\n", bbox[0], bbox[3]);
		printf("bbox: y:  %f - %f\n", bbox[1], bbox[4]);
		printf("bbox: z:  %f - %f\n", bbox[2], bbox[5]);
        float miny = bbox[1];
        float maxy = bbox[4];
        
        float width = atof(argv[5]);
        if (argc>6)  { mold_mode = atoi(argv[6]);  }
        MeshTriangles* mt = parse_triangles(argv[2],(float) width);
        if (mt == NULL) {return 1;} 
        printf("parsed interpolation\n");

        printf("n points=%d\n",mt->npoints);
        printf("n triangles=%d\n",mt->ntriangles);
        meshindex* mi = build_mesh_index( mt->points, mt->triangles, mt->ntriangles, 100, 100 ) ;


        MeshTriangles* texture = parse_triangles_with_normals(argv[1],(float) 1.0);
        if (texture == NULL) {return 1;} 

/*
        for(int i = 0;i<texture->ntriangles;i++) {
            for (int j=0;j<3;j++) {
                int a = texture->triangles[i*3 + j];
                float x = texture->xpoints[a*3];
                float y = texture->xpoints[a*3+1];
                float z = texture->xpoints[a*3+2];

                meshindex_it* it = index_iterator( x,y,mi );
                printf("\n");
                while (!it->done) { 
                    int aa = mt->triangles[3*it->t    ];
                    int bb = mt->triangles[3*it->t + 1];
                    int cc = mt->triangles[3*it->t + 2];
                    printf("triangle %d (%f,%f) index %d   (%f,%f)-(%f,%f)-(%f,%f)\n",i,x,y,it->t, 
                        mt->points[2*aa], 
                        mt->points[2*aa+1], 
                        mt->points[2*bb], 
                        mt->points[2*bb+1], 
                        mt->points[2*cc], 
                        mt->points[2*cc+1]
                    );
                    it = next_index(it);
                }
            }
            
        }
*/

        printf("parsed with normals\n");  fflush(stdout);
        apply_interpolation_to_mesh( texture, mt, 10.0 );

        free_meshtriangles( mt ) ;

        //apply_transform_to_mesh( texture, transf_X );
        //write_to_stl(texture,stlfile);

        if (! mold_mode) {
            float thickness = 15.0;



            float offset[3] = { 0,0,0};
            Plane clip_rs_planes[3];
            clip_rs_planes[0] = (Plane){ {0,  0, 1},  0.1  };

            Plane clip_xyz_planes[2];
            clip_xyz_planes[0] = (Plane){ {0,1,0}, 0.0 };
            clip_xyz_planes[1] = (Plane){ {0,-1,0}, bbox[4] };
            

            for (int i=0; i<texture->ntriangles; i++) {
                //printf("mt path %d\n",i);
                
                int n1 = texture->triangles[3*i];
                int n2 = texture->triangles[3*i+1];
                int n3 = texture->triangles[3*i+2];

                if ( texture->nxpoints <= n1 ) {printf("n1=%d but  texture->nxpoints=%d i=%d ntriangles=%d\n",n1, 
                texture->nxpoints,i, texture->ntriangles); fflush(stdout);;}

                if ( texture->nxpoints <= n2 ) {printf("n1=%d but  texture->nxpoints=%d i=%d ntriangles=%d\n",n2, 
                texture->nxpoints,i, texture->ntriangles); fflush(stdout);;}

                if ( texture->nxpoints <= n3 ) {printf("n1=%d but  texture->nxpoints=%d i=%d ntriangles=%d\n",n3, 
                texture->nxpoints,i, texture->ntriangles); fflush(stdout);}

                assert( texture->nxpoints > n1 );
                assert( texture->nxpoints > n2 );
                assert( texture->nxpoints > n3 );

                //tessAddContour(tess, 2, &path_points[ path_offsets[i] ] , sizeof(float)*2, path_lengths[i]);
                float v1[3] = { texture->xpoints[n1*3+0] ,texture->xpoints[n1*3+1] ,texture->xpoints[n1*3+2] };
                float v3[3] = { texture->xpoints[n2*3+0] ,texture->xpoints[n2*3+1] ,texture->xpoints[n2*3+2]};
                float v2[3] = { texture->xpoints[n3*3+0] ,texture->xpoints[n3*3+1] ,texture->xpoints[n3*3+2] };
    
    //            assert( texture->nxpoints < n1 );
    //            float v2[3] = { texture->paths[i*6+2] ,texture->paths[i*6+3] , 0.1 };
    //            float v3[3] = { texture->paths[i*6+4] ,texture->paths[i*6+5] , 0.1 };

                if (v1[2]<0.0 && v2[2]<0.0 && v3[2]<0.0) {continue;}

                printf("print clipped %d %d %d\n",n1,n2,n3);

                print_triangle_clipped(v1,v2,v3, stlfile, clip_rs_planes, 1, clip_xyz_planes, 2, &offset[0],0);
                            fflush( stlfile );
                if (i%1000==0) {printf("%d/%d\n",i,texture->ntriangles); fflush(stdout);}
      //          print_triangle_clipped(v1,v2,v3, stlfile, clip_rs_planes, 2, clip_xyz_planes, 2,&offsets[quadrant*3]); //2 omits inside clip
    //            print_triangle_bbox(v1,v2,v3, stlfile, bbox);
            }
            printf("done writing texture\n");
            fflush(stdout);

//            write_surface_stl( &transform, stlfile, 0, 1, thickness,&offset[0] );
//            write_texture_back_stl( &transform, stlfile, 0, 1, 0.5*thickness,&offset[0]) ;


            fprintf(stlfile, "endsolid x\n");
            fclose(stlfile);
            printf("done closing stlfile\n"); fflush(stdout);
            free_meshtriangles( texture ) ;

            return(0);
        } 

            float offset_size = 50.0;
            float thickness = 20.0;
            float angle_pad = 0.005;


        float offsets[18] = {
            offset_size,0,0,
            0,0,offset_size,
            -offset_size,0,0,
            0,0,-offset_size,
            offset_size,0,0,
            0,0,offset_size
        };

//        for(int quadrant=1; quadrant<2; quadrant++ ) {
        for(int quadrant=0; quadrant<6; quadrant++ ) {

            float offset[3];  
            float mins = 0; //= atan( spec.squash)/(2.0*M_PI) ;
            float maxs = 0; // = 0.5- (atan( spec.squash ) / (2.0*M_PI));


            set(offset, &offsets[quadrant*3] );

            switch (quadrant) {
                case 0:
                    //set(offset, (float[3]) { offset_size,0,0 });
                    mins = 0.0 - (atan( spec.squash ) / (2.0*M_PI)) - angle_pad;
                    maxs = 0.0 + (atan( spec.squash)/(2.0*M_PI)) + angle_pad ;

                    break;
                case 1:
                    //set(offset, (float[3]) { 0, 0,offset_size });
                    mins = atan( spec.squash)/(2.0*M_PI) - angle_pad;
                    maxs = 0.5 - (atan( spec.squash ) / (2.0*M_PI)) + angle_pad;

                    break;
                case 2:
                    //set(offset, (float[3]) { -offset_size,0,0 });
                    mins = 0.5 - (atan( spec.squash ) / (2.0*M_PI))  - angle_pad; // atan( spec.squash)/(2.0*M_PI) ;
                    maxs = 0.5 + (atan( spec.squash ) / (2.0*M_PI)) + angle_pad;

                    break;
                case 3:
                    //set(offset, (float[3]) { 0,0,-offset_size });
                    mins = 0.5 + (atan( spec.squash ) / (2.0*M_PI)) - angle_pad; // atan( spec.squash)/(2.0*M_PI) ;
                    maxs = 1.0 - (atan( spec.squash ) / (2.0*M_PI)) + angle_pad;

                    break;
                case 4:
                    //set(offset, (float[3]) { offset_size,0,0 });
                    mins = 1.0 - (atan( spec.squash ) / (2.0*M_PI)) - angle_pad;
                    maxs = 1.0 + (atan( spec.squash)/(2.0*M_PI))  + angle_pad;

                    break;
                case 5:
                    //set(offset, (float[3]) { offset_size,0,0 });
                    mins = 1+atan( spec.squash)/(2.0*M_PI) - angle_pad;
                    maxs = 1+0.5 - (atan( spec.squash ) / (2.0*M_PI)) + angle_pad;

                    break;
            }

            printf("quadrant %d:  offset %f %f %.0f %.0f %.0f\n",quadrant,mins,maxs, offset[0],offset[1],offset[2]);

            printf("atan thing = %f\n", -mins  );
            printf("atan thing = %f\n", maxs  );

            Plane clip_rs_planes[3];
            clip_rs_planes[0] = (Plane){ {0,  1, 0}, -mins };
            clip_rs_planes[1] = (Plane){ {0, -1, 0},  maxs  };
            clip_rs_planes[2] = (Plane){ {0,  0, 1},  0.1  };

            Plane clip_xyz_planes[2];
            clip_xyz_planes[0] = (Plane){ {0,1,0}, 0.0 };
            clip_xyz_planes[1] = (Plane){ {0,-1,0}, bbox[4] };
            

            for (int i=0; i<texture->ntriangles; i++) {
                //printf("mt path %d\n",i);
                int n1 = texture->triangles[3*i];
                int n2 = texture->triangles[3*i+1];
                int n3 = texture->triangles[3*i+2];

                if ( texture->nxpoints <= n1 ) {printf("n1=%d but  texture->nxpoints=%d i=%d ntriangles=%d\n",n1, 
                texture->nxpoints,i, texture->ntriangles);}
                assert( texture->nxpoints > n1 );
                assert( texture->nxpoints > n2 );
                assert( texture->nxpoints > n3 );

                //tessAddContour(tess, 2, &path_points[ path_offsets[i] ] , sizeof(float)*2, path_lengths[i]);
                float v1[3] = { texture->xpoints[n1*3+0] ,texture->xpoints[n1*3+1] ,texture->xpoints[n1*3+2] };
                float v3[3] = { texture->xpoints[n2*3+0] ,texture->xpoints[n2*3+1] ,texture->xpoints[n2*3+2]};
                float v2[3] = { texture->xpoints[n3*3+0] ,texture->xpoints[n3*3+1] ,texture->xpoints[n3*3+2] };
    
    //            assert( texture->nxpoints < n1 );
    //            float v2[3] = { texture->paths[i*6+2] ,texture->paths[i*6+3] , 0.1 };
    //            float v3[3] = { texture->paths[i*6+4] ,texture->paths[i*6+5] , 0.1 };

                if (v1[2]<0.0 && v2[2]<0.0 && v3[2]<0.0) {continue;}

                print_triangle_clipped(v1,v2,v3, stlfile, clip_rs_planes, 3, clip_xyz_planes, 2,offset, -0.5*thickness);
                            fflush( stlfile );

      //          print_triangle_clipped(v1,v2,v3, stlfile, clip_rs_planes, 2, clip_xyz_planes, 2,&offsets[quadrant*3]); //2 omits inside clip
    //            print_triangle_bbox(v1,v2,v3, stlfile, bbox);
            }

            fflush( stlfile );
            printf("done writing triangles for texture\n");
            fflush(stdout);

            if (quadrant>=4) {continue;}

            float length = 300.0;
            write_surface_stl( &transform, stlfile, mins, maxs, thickness,&offsets[quadrant*3] );
            write_texture_back_stl( &transform, stlfile, mins, maxs, 0.5*thickness,&offsets[quadrant*3]) ;

          //  if (quadrant==1) {
             //   write_floor_flange_stl( &spec, stlfile, mins, maxs,length, thickness/(bbox[4]-bbox[2]),offset, 0 ,0);
                
                // too overhanging for 3d print?
//                write_floor_flange_stl( &spec, stlfile, mins, maxs,length, thickness/(bbox[4]-bbox[2]),offset, 1 ,1);
         //   }
 

            stl_output_config outputcf;
            outputcf.fp = stlfile;
            outputcf.clipping_planes = &clip_xyz_planes[0];
            outputcf.n_clipping_planes = 2;

            write_parting_sufrace_stl( quadrant,   length , &spec,  outputcf, &offsets[quadrant*3], thickness, 0); 
            write_parting_sufrace_stl( quadrant-1, length , &spec,  outputcf, &offsets[quadrant*3], thickness, 1); 

            write_support_ties_stl(&spec, quadrant, quadrant, 
            &offsets[3* quadrant ], 
            &offsets[3* ((quadrant+1)%4)], length, thickness, outputcf );
    
        }

	fprintf(stlfile, "endsolid x\n");
    fclose(stlfile);
        free_meshtriangles( texture ) ;

    }

    return(0);
}

