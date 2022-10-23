#include <stdlib.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "surface.h"

#include "frozen.h"
#include "catmullrom.h"
#include "vectorops.h"
#include "triangles.h"
#include "assert.h"



void wrinkley(float* xyz, float* vtheta) {

    double v = (double) vtheta[0];
    double theta = (double) vtheta[1];

    double phi = M_PI*(0.6+0.6*(v-0.5));

    double h0 = 1.7;

    double h = -h0*cos(phi);
    double r= (1.1-0.1*pow(sin(theta*7.0 - M_PI*v*2.5),6.0) )*sin(phi);

//    double r1 = 1.5;
//    double r2 = 9.0;
    //double ee=1.3;
    double ee=1.05;

    //double eeh=2.0;

//    double r = r1 + v*(r2-r1);

    xyz[0] = ee*r*sin( theta  );  // X
    xyz[1] = h;  // Y
    xyz[2] = r*cos( theta  );  // Z
}

void flattened_vase(float* xyz, float* vtheta) {
    double v = (double) vtheta[0];
    double theta = (double) vtheta[1];
    double h0 = 10.0;
    double r1 = 1.5;
    double r2 = 9.0;
    double ee=1.6;

    double r = r1 + v*(r2-r1);

    xyz[0] = (float) sqrt(r)*sin( theta  );  // X
    xyz[1] = (float) v*h0;  // Y
    xyz[2] = (float) ee*sqrt(r)*cos( theta );  // Z

   // printf("X: %f %f -> %f %f %f\n",vtheta[0], vtheta[1], xyz[0], xyz[1], xyz[2]);

}

void bulbous_vase(float* xyz, float* vtheta) {
    double v = (double) vtheta[0];
    double theta = (double) vtheta[1];
    double h0 = 10.0;
    double r0=6.0;
    double r1 = 1.5;
    double r2 = 9.0;
    double ee=1.2;

    double phi = M_PI*(0.6+0.6*(0.5-v));
    double rxy =r0*sin(phi);// r1 + v*(r2-r1);

    xyz[0] = (float) rxy*sin( theta  );  // X
    xyz[1] = (float) r0*cos(phi);  // Y
    xyz[2] = (float) ee*rxy*cos( theta );  // Z

   // printf("X: %f %f -> %f %f %f\n",vtheta[0], vtheta[1], xyz[0], xyz[1], xyz[2]);

}

void bulbous_vase2(float* xyz, float* vtheta) {
    double v = (double) vtheta[0];
    double theta = (double) vtheta[1];
    double h0 = 10.0;
    double r0=6.0;
    double r1 = 1.5;
    double r2 = 9.0;
    double ee=1.2;

    double phi = M_PI*(0.6+0.6*(0.5-v));
    double rxy =r0*sin(phi);// r1 + v*(r2-r1);

    xyz[0] = (float) sqrt(rxy)*sin( theta  );  // X
    xyz[1] = (float) r0*cos(phi);  // Y
    xyz[2] = (float) ee*sqrt(rxy)*cos( theta );  // Z

   // printf("X: %f %f -> %f %f %f\n",vtheta[0], vtheta[1], xyz[0], xyz[1], xyz[2]);

}



void X(float* xyz, float* vtheta) {
    //flattened_vase(xyz, vtheta);
    bulbous_vase2(xyz, vtheta);
    //wrinkley(xyz, vtheta);
    
}

//void Xs(float* xyz, float* uv, pottoy_spec_t* spec ) {
//}
/*void spline_getPoint( float t, float* point,pottoy_spec_t*  spec ) {
}*/

void sor(float* xyz, float* uv, pottoy_spec_t* spec ) {

    float point[2];
//	spline_getPoint( uv[0], &point[0], spec );

    if (uv[0]>1.0) {uv[0]=1.0;}
    if (uv[0]<0.0) {uv[0]=0.0;}

    catmullrom2(uv[0],&spec->points[0], spec->npoints, &point[0]);

	float r = point[0];
	float y = point[1];

    xyz[0] = r*cos( uv[1]* M_PI *2.0) ;
    xyz[1] = y;
    xyz[2] = spec->squash*r*sin( uv[1] * M_PI *2.0);

}


#define TBUFF_SIZE2 3*3*32*2
void output_triangle( float* v1, float* v2, float* v3, stl_output_config cf ) {

    float triangles[ TBUFF_SIZE2 ];
  
    int nt=1;
    int ti=0;
    int r;
    int altered;
    int nt2=1;
    int ti2=0;

//    float drop_offset[3] = { 0,0,-dropdepth };

    set( &triangles[0], v1 );
    set( &triangles[3], v2 );
    set( &triangles[6], v3 );
  //  printf("clip triangle\n "); fflush(stdout);

    int n_written=0;
    while(ti2<nt2) {
        int nr = 1;
        int nt20 = nt2;
        for(int pi=0;pi<cf.n_clipping_planes ;pi++) {
            //printf("clip with xyzplane ");
           // printf_triangle(&triangles2[ ti2*9 ]);
            int tr = 0; 
            clip_triangle( &triangles[ ti2*9 ] , &triangles[ ti2*9 +3], &triangles[ ti2*9 +6], 
            &triangles[ nt2*9 ], &cf.clipping_planes[pi], &tr, &altered );
            nr *= tr;
            if (tr>1) { 
              //  printf_triangle(&triangles2[ nt2*9 ]);
            nt2++; //printf("nt=%d\n",nt);
             assert(nt*9<TBUFF_SIZE2);} 
        }
        if (nr>0) {
               
            stl_triangle( 
                &triangles[ ti2*9 ]  , 
                &triangles[ ti2*9 +3] , 
                &triangles[ ti2*9 +6] , cf.fp );
               // fflush(fp);

            n_written++;
        } 

        ti2++;
    }
 //   stl_triangle( v1, v2, v3, cf.fp );
}


void stl_triangle( float* v1, float* v2, float* v3, FILE* fp ) {
    float n[3];

    triangle_normal(v1,v2,v3,&n[0]);

	fprintf(fp,"  facet normal %f %f %f\n", n[0], n[1], n[2]);
	fprintf(fp,"    outer loop\n");
	fprintf(fp,"      vertex %f %f %f\n",v1[0],-v1[2],v1[1]);
	fprintf(fp,"      vertex %f %f %f\n",v2[0],-v2[2],v2[1]);
	fprintf(fp,"      vertex %f %f %f\n",v3[0],-v3[2],v3[1]);

	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");

}


void output_quad( float* ray1 , float* ray2 ,  stl_output_config cf , int flip ) {

    if (flip) {
        output_triangle( &ray1[0], &ray1[3], &ray2[3] ,cf );
        output_triangle( &ray1[0], &ray2[3], &ray2[0] ,cf );
    } else {
        output_triangle( &ray1[0], &ray1[3], &ray2[0] ,cf );
        output_triangle( &ray2[0], &ray1[3], &ray2[3] ,cf );
    }

}

void stl_printquad( float* ray1 , float* ray2 , FILE* fp, int flip ) {

    if (flip) {
        stl_triangle( &ray1[0], &ray1[3], &ray2[3] ,fp );
        stl_triangle( &ray1[0], &ray2[3], &ray2[0] ,fp );
    } else {
        stl_triangle( &ray1[0], &ray1[3], &ray2[0] ,fp );
        stl_triangle( &ray2[0], &ray1[3], &ray2[3] ,fp );
    }
}

void mold_parting_surface_norm(float u, int quadrant, pottoy_spec_t* spec, float* result) ;




void output_cuboid(float* c, stl_output_config cf) {

    output_triangle( &c[0],  &c[6] , &c[3] ,cf );  //edge of the thing
    output_triangle( &c[6],  &c[9] , &c[3] ,cf );  //edge of the thing

    output_triangle( &c[12],  &c[15] , &c[18] ,cf );  //edge of the thing
    output_triangle( &c[15],  &c[21] , &c[18] ,cf );  //edge of the thing

    output_triangle( &c[0],  &c[12] , &c[6] ,cf );  //edge of the thing
    output_triangle( &c[6],  &c[12] , &c[18] ,cf );  //edge of the thing

    output_triangle( &c[3],  &c[9] , &c[15] ,cf );  //edge of the thing
    output_triangle( &c[15],  &c[9] , &c[21] ,cf );  //edge of the thing

    output_triangle( &c[3],  &c[15] , &c[0] ,cf );  //edge of the thing
    output_triangle( &c[0],  &c[15] , &c[12] ,cf );  //edge of the thing

    output_triangle( &c[6],  &c[18] , &c[9] ,cf );  //edge of the thing
    output_triangle( &c[9],  &c[18] , &c[21] ,cf );  //edge of the thing

}

void stl_cuboid(float* c, FILE* stlfile) {

    stl_triangle( &c[0],  &c[6] , &c[3] ,stlfile );  //edge of the thing
    stl_triangle( &c[6],  &c[9] , &c[3] ,stlfile );  //edge of the thing

    stl_triangle( &c[12],  &c[15] , &c[18] ,stlfile );  //edge of the thing
    stl_triangle( &c[15],  &c[21] , &c[18] ,stlfile );  //edge of the thing

    stl_triangle( &c[0],  &c[12] , &c[6] ,stlfile );  //edge of the thing
    stl_triangle( &c[6],  &c[12] , &c[18] ,stlfile );  //edge of the thing

    stl_triangle( &c[3],  &c[9] , &c[15] ,stlfile );  //edge of the thing
    stl_triangle( &c[15],  &c[9] , &c[21] ,stlfile );  //edge of the thing

    stl_triangle( &c[3],  &c[15] , &c[0] ,stlfile );  //edge of the thing
    stl_triangle( &c[0],  &c[15] , &c[12] ,stlfile );  //edge of the thing

    stl_triangle( &c[6],  &c[18] , &c[9] ,stlfile );  //edge of the thing
    stl_triangle( &c[9],  &c[18] , &c[21] ,stlfile );  //edge of the thing

}


float hole_width(float t, int n_holes, float size) {

    float tp = t*n_holes; 
    int n = (int) tp;
    tp = (tp - (float) n)-0.5;

    if (fabs(tp)>0.2) {return(0.0);}
    return(0.5*size);

}


void write_support_ties_stl(
    pottoy_spec_t* spec,
    int quadrant1, int quadrant2, float* offset1, 
        float* offset2, float length, float thickness, stl_output_config cf ) {

    int n_levels = 20;
    int n_per_level = 4;

    float th = 0.025;
    float wi = 0.0025;

    float ray1[6];
    float ray2[6];

    float ray3[6];
    float ray4[6];

    float cuboid[24];

    float delta_t = 1.0/((float) n_levels);

    for (int j=0; j<n_per_level; j++) { 
        float r_fact = ((float) j) / ((float)(n_per_level + 1));
    for( int i = 0; i<n_levels; i++ ) {
        float t = i*delta_t;

        if ( hole_width(t,3,0.3*length)>0 ) {continue;} 

        mold_parting_norm_ray( t, length , thickness, quadrant1, spec, &ray1[0]);
        add(ray1,ray1,offset1);
        add(&ray1[3],&ray1[3],offset1);

        mold_parting_norm_ray( t, length , thickness, quadrant2, spec, &ray2[0]);
        add(ray2,     ray2,   offset2);
        add(&ray2[3],&ray2[3],offset2);

        mold_parting_norm_ray( t+wi, length , thickness, quadrant1, spec, &ray3[0]);
        add(ray3,ray3,offset1);
        add(&ray3[3],&ray3[3],offset1);

        mold_parting_norm_ray( t+wi, length , thickness, quadrant2, spec, &ray4[0]);
        add(ray4,ray4,offset2);
        add(&ray4[3],&ray4[3],offset2);

//        set(           &cuboid[0], &ray1[3] );
        weighted_sum3( &cuboid[0], (1-r_fact), r_fact, &ray1[3], &ray1[0] );

//        set(           &cuboid[3], &ray3[3] );
        weighted_sum3( &cuboid[3], (1-r_fact-th), r_fact+th, &ray1[3], &ray1[0] );

        //set(           &cuboid[6], &ray2[3] );
        weighted_sum3( &cuboid[6], (1-r_fact), r_fact, &ray2[3], &ray2[0] );

        weighted_sum3( &cuboid[9], (1-r_fact-th), r_fact+th, &ray2[3], &ray2[0] );

        //set(           &cuboid[12], &ray3[3] );
        weighted_sum3( &cuboid[12], (1-r_fact), r_fact, &ray3[3], &ray3[0] );
        weighted_sum3( &cuboid[15], (1-r_fact-th), r_fact+th, &ray3[3], &ray3[0] );
        //set(           &cuboid[18], &ray4[3] );
        weighted_sum3( &cuboid[18], (1-r_fact), r_fact, &ray4[3], &ray4[0] );
        weighted_sum3( &cuboid[21], (1-r_fact-th), r_fact+th, &ray4[3], &ray4[0] );


        output_cuboid( cuboid, cf );
//        stl_triangle( &ray1[3],  &ray2[0] , &ray2[0] ,stlfile );  //edge of the thing

        //mold_parting_norm_ray( t, length , thickness, quadrant, spec, &ray1[0]);
    }}

}

/*
void write_support_ties_stl(int quadrant, float length, float* offset, float thickness,FILE* stlfile ) {

    int n_levels = 10;
    int n_per_level = 4;

    float ray1[6];
    float delta_t = 1.0/((float) n_levels);

    for( int i = 0; i<n_levels; i++ ) {
        float t = i*delta_t;

        //mold_parting_norm_ray( t, length , thickness, quadrant, spec, &ray1[0]);
    }

}
*/


void output_split_quad( float* ray1 , float* ray2 ,stl_output_config cf , int flip, float a, float b, float c, float d ) {

    float ray1a[6]; 
    float ray1b[6]; 

    float ray2a[6]; 
    float ray2b[6]; 

    weighted_sum3( ray1a   , 1  ,  0 , ray1 , &ray1[3]);	
    weighted_sum3(&ray1a[3], 1-a,  a , ray1 , &ray1[3]);	

    weighted_sum3( ray1b   , 1-b ,  b , ray1 , &ray1[3]);	
    weighted_sum3(&ray1b[3], 0   ,  1 , ray1 , &ray1[3]);	

    weighted_sum3( ray2a   , 1  ,  0 , ray2 , &ray2[3]);	
    weighted_sum3(&ray2a[3], 1-c,  c , ray2 , &ray2[3]);	

    weighted_sum3( ray2b   , 1-d ,  d , ray2 , &ray2[3]);	
    weighted_sum3(&ray2b[3], 0   ,  1 , ray2 , &ray2[3]);	

    output_quad( ray1a , ray2a , cf, flip ) ;
    output_quad( ray1b , ray2b , cf, flip ) ;

}


void stl_print_split_quad( float* ray1 , float* ray2 , FILE* fp, int flip, float a, float b, float c, float d ) {

    float ray1a[6]; 
    float ray1b[6]; 

    float ray2a[6]; 
    float ray2b[6]; 

    weighted_sum3( ray1a   , 1  ,  0 , ray1 , &ray1[3]);	
    weighted_sum3(&ray1a[3], 1-a,  a , ray1 , &ray1[3]);	

    weighted_sum3( ray1b   , 1-b ,  b , ray1 , &ray1[3]);	
    weighted_sum3(&ray1b[3], 0   ,  1 , ray1 , &ray1[3]);	

    weighted_sum3( ray2a   , 1  ,  0 , ray2 , &ray2[3]);	
    weighted_sum3(&ray2a[3], 1-c,  c , ray2 , &ray2[3]);	

    weighted_sum3( ray2b   , 1-d ,  d , ray2 , &ray2[3]);	
    weighted_sum3(&ray2b[3], 0   ,  1 , ray2 , &ray2[3]);	

    stl_printquad( ray1a , ray2a , fp, flip ) ;
    stl_printquad( ray1b , ray2b , fp, flip ) ;

}


void write_parting_sufrace_stl( int quadrant , float l, pottoy_spec_t* spec, stl_output_config cf, float* offset, 
      float thickness, int reverse) {

   // FILE* stlfile = cf.fp;

    float ray1[6];
    float ray2[6];

    float ray1_back[6];
    float ray2_back[6];
    float norm[3];

    float hole_contour_points1[256];
    float hole_contour_points2[256];
    int hci1=0;
    int hci2=0;
    int in_hole=0;

    int n_segments = 51;
    float delta_t = 1.0 / ((float) n_segments);
    float t=0.0;

    float gap1, gap2;
    int n_keys = 3;

    mold_parting_norm_ray( t, l , thickness, quadrant, spec, &ray1[0]);
    add(ray1,ray1,offset);
    add(&ray1[3],&ray1[3],offset);

    mold_parting_surface_norm(t, quadrant,  spec, norm); scale(norm, reverse ? thickness : -thickness);
    memcpy( &ray1_back[0] , &ray1[0] , sizeof(float)*6 );
    add(ray1_back,ray1_back, norm );
    add(&ray1_back[3],&ray1_back[3], norm );
    gap1 = hole_width(t,n_keys,0.3*l);
    if (gap1>0) { in_hole=1; }

    if (ray1[1] < ray1[4] ) { 
        float pv1[3];
        float pvf[3];
        float pvb[3];
        diff3(&ray1[3],ray1,pv1);
        pv1[1]=0.0;
        normalize(pv1,pv1);
        scale(pv1,l);
        add(pvf,pv1,ray1);
        add(pvb,pv1,ray1_back);

//        float myhat[3] = { 0,-l,0 };
//        add(pv1,ray1,myhat);
//        add(pvb,ray1_back,myhat);
        printf("wedge needed at bottom: %f < %f\n", ray1[1], ray1[4]); 

        if (reverse) {
            output_triangle( &ray1[0],pvf,&ray1[3] ,cf ); 
            output_triangle( &ray1_back[0],&ray1_back[3],pvb ,cf ); 

            output_triangle( pvf , pvb, &ray1[3],cf ); 
            output_triangle( &ray1[3],pvb,&ray1_back[3],cf ); 
        } else {
            output_triangle( &ray1[0],&ray1[3],pvf ,cf ); 
            output_triangle( &ray1_back[0],pvb,&ray1_back[3] ,cf ); 

            output_triangle( pvb , pvf, &ray1[3],cf ); 
            output_triangle( pvb,&ray1[3],&ray1_back[3],cf ); 
//            output_triangle( &ray1_back[0],pvb,&ray1_back[3] ,cf ); 


        }
    }

    for(int i=0; i<n_segments; i++) {
       // printf("parting surface i=%d t=%f ( %f %f %f => %f %f %f\n",i,t, 
        //    ray1[0], ray1[1], ray1[2],  ray1[3], ray1[4], ray1[5] );

        t += delta_t;
        gap2 = hole_width(t,n_keys,0.3*l);

        mold_parting_norm_ray( t, l, thickness , quadrant, spec, &ray2[0]);
        add(ray2,ray2,offset);
        add(&ray2[3],&ray2[3],offset);

        mold_parting_surface_norm(t, quadrant,  spec, norm); scale(norm,reverse ? thickness : -thickness);
        memcpy( &ray2_back[0] , &ray2[0] , sizeof(float)*6 );
        add(ray2_back, ray2_back, norm );
        add(&ray2_back[3], &ray2_back[3], norm );

        if (reverse) {
            if ((gap1==0) && (gap2 ==0)) {
                output_quad( &ray1[0],&ray2[0],cf, 0 );
                output_quad( &ray2_back[0],&ray1_back[0],cf, 1 );
            } else {
                output_split_quad(&ray1[0],&ray2[0],cf, 0 , 0.5-gap1/l, 0.5+gap1/l, 0.5-gap2/l, 0.5+gap2/l); 
                output_split_quad(&ray2_back[0],&ray1_back[0],cf, 1 , 0.5-gap2/l, 0.5+gap2/l,0.5-gap1/l, 0.5+gap1/l); 
            }

            output_triangle( &ray1[3], &ray1_back[3], &ray2[3] ,cf );  //edge of the thing
            output_triangle( &ray2[3], &ray1_back[3], &ray2_back[3],cf );  //edge of the thing
            
            output_triangle( &ray1[0],  &ray2[0] ,&ray1_back[0],cf );  //edge of the thing
            output_triangle( &ray2[0], &ray2_back[0], &ray1_back[0],cf );  //edge of the thing
//            stl_triangle( &ray2[0], &ray1[3], &ray2[3] ,fp ); 

        } else {
            if ((gap1==0) && (gap2 ==0)) {
                output_quad( &ray2[0],&ray1[0],cf, 1 );
                output_quad( &ray1_back[0],&ray2_back[0],cf, 0 );
            } else {
                output_split_quad(&ray2[0],&ray1[0],cf, 0,  0.5-gap2/l, 0.5+gap2/l,0.5-gap1/l, 0.5+gap1/l ); 
                output_split_quad(&ray1_back[0],&ray2_back[0],cf, 1 , 0.5-gap1/l, 0.5+gap1/l, 0.5-gap2/l, 0.5+gap2/l ); 
            }

            output_triangle( &ray1[3], &ray2[3], &ray1_back[3] ,cf );  //edge of the thing
            output_triangle( &ray2[3], &ray2_back[3], &ray1_back[3],cf );  //edge of the thing

            output_triangle( &ray1[0], &ray1_back[0], &ray2[0] ,cf );  //edge of the thing
            output_triangle( &ray2[0], &ray1_back[0], &ray2_back[0],cf );  //edge of the thing

          //  stl_triangle( &ray1[0], &ray1[3], &ray2[0] ,fp );  //edge of the thing
          //  stl_triangle( &ray2[0], &ray1[3], &ray2[3] ,fp ); 

        }

        if ((gap1==0) && (gap2==0)) {
            if (hci1>0 || hci2>0) {

                float key_triangle_buffer[300];
                float key_offset[3];

                set(key_offset,norm);

                if ((quadrant+reverse)%2==0) scale(key_offset,-1.0);
                // make the key here

                
                while (hci2 >0 ) {
                    set( &hole_contour_points1[(hci1++)*3], &hole_contour_points2[(--hci2)*3] );
                }

                int nt_key = contour_to_mesa(hci1, &hole_contour_points1[0], key_offset, 0.5, 100, &key_triangle_buffer[0]);

                if (reverse) {
                    for (int i=0;i<nt_key;i++) {
                        output_triangle( &key_triangle_buffer[9*i],&key_triangle_buffer[9*i+6],&key_triangle_buffer[9*i+3],cf);   
                    }
                } else {
                    for (int i=0;i<nt_key;i++) {
                        output_triangle( &key_triangle_buffer[9*i],&key_triangle_buffer[9*i+3],&key_triangle_buffer[9*i+6],cf);   
                    }
                }

                for (int i=0;i<nt_key*3;i++) {
                    add( &key_triangle_buffer[3*i],&key_triangle_buffer[3*i],norm );
                }

                if (!reverse) {
                    for (int i=0;i<nt_key;i++) {
                        output_triangle( &key_triangle_buffer[9*i],&key_triangle_buffer[9*i+6],&key_triangle_buffer[9*i+3],cf);   
                    }
                } else {
                    for (int i=0;i<nt_key;i++) {
                        output_triangle( &key_triangle_buffer[9*i],&key_triangle_buffer[9*i+3],&key_triangle_buffer[9*i+6],cf);   
                    }
                }

                hci1=0;
                hci2=0;
                in_hole=0;
            }
        } else {
//            if (in_hole==0) {
            if ((gap1==0) && (gap2>0)  ) {
                // start a new hole


                weighted_sum3(  &hole_contour_points1[(hci1++)*3]   , 0.5  ,  0.5 , ray1 , &ray1[3]);	 // apex point

                weighted_sum3(  &hole_contour_points1[(hci1++)*3]   , 1-(0.5-gap2/l)  ,  (0.5-gap2/l) , ray2 , &ray2[3]);	 // apex point
                weighted_sum3(  &hole_contour_points2[(hci2++)*3]   , 1-(0.5+gap2/l)  ,  (0.5+gap2/l), ray2 , &ray2[3]);	 // apex point


                in_hole = 1;
            } else if (gap2==0.0) {   // end a hole

                weighted_sum3(  &hole_contour_points1[(hci1++)*3]   , 0.5  ,  0.5 , ray2 , &ray2[3]);	 // apex point
                

            } else  {  // continue a hole
                weighted_sum3(  &hole_contour_points1[(hci1++)*3]   , 1-(0.5-gap2/l)  ,  (0.5-gap2/l) , ray2 , &ray2[3]);	 // apex point
                weighted_sum3(  &hole_contour_points2[(hci2++)*3]   , 1-(0.5+gap2/l)  ,  (0.5+gap2/l), ray2 , &ray2[3]);	 // apex point

            }

            

        }

//        if (i==n_segments-1) {        }

        memcpy( &ray1[0] , &ray2[0] , sizeof(float)*6 );
        memcpy( &ray1_back[0] , &ray2_back[0] , sizeof(float)*6 );
        gap1=gap2;
    }

    if (ray1[1] > ray1[4] ) { 
        float pv1[3];
        float pvf[3];
        float pvb[3];
        diff3(&ray1[3],ray1,pv1);
        pv1[1]=0.0;
        normalize(pv1,pv1);
        scale(pv1,l);
        add(pvf,pv1,ray1);
        add(pvb,pv1,ray1_back);

//        float myhat[3] = { 0,-l,0 };
//        add(pv1,ray1,myhat);
//        add(pvb,ray1_back,myhat);
        printf("wedge needed at bottom: %f < %f\n", ray1[1], ray1[4]); 

        if (!reverse) {
            output_triangle( &ray1[0],pvf,&ray1[3] ,cf ); 
            output_triangle( &ray1_back[0],&ray1_back[3],pvb ,cf ); 

            output_triangle( pvf , pvb, &ray1[3],cf ); 
            output_triangle( &ray1[3],pvb,&ray1_back[3],cf ); 
        } else {
            output_triangle( &ray1[0],&ray1[3],pvf ,cf ); 
            output_triangle( &ray1_back[0],pvb,&ray1_back[3] ,cf ); 

            output_triangle( pvb , pvf, &ray1[3],cf ); 
            output_triangle( pvb,&ray1[3],&ray1_back[3],cf ); 
//            output_triangle( &ray1_back[0],pvb,&ray1_back[3] ,cf ); 


        }
    }
    


}


void xz_project_unit( float* xyz ) {
    xyz[1]=0.0;
    normalize(xyz,xyz);
}

void unfaceted_surface_norm( float t, float s, pottoy_spec_t* spec, float* result ) {

  //  float curve_p[2];
  //  catmullrom2(u,&spec->points[0], spec->npoints, &curve_p[0]);


    double theta = s* 2.0*M_PI ;
    double cos_theta=cos(theta);
    double sin_theta=sin(theta);

    float t1[3];  // 45 degrees to x and y
    float t2rs[2];  // along the spline
    float t2[3];  // along the spline



    t1[0] = -sin_theta  ;
    t1[1] =  0.0;
    t1[2] =  spec->squash*cos_theta;

    catmullrom2_tangent(t, &spec->points[0], spec->npoints, &t2rs[0]);

    t2[0] = t2rs[0] * cos_theta ;
    t2[1] = t2rs[1];
    t2[2] = t2rs[0] * spec->squash*sin_theta;

    float nn[3];
//    float n[3];

    cross_product(&t2[0],&t1[0],&nn[0]);
    normalize(&nn[0],result);

}


void write_floor_flange_stl(  pottoy_spec_t* spec, FILE* stlfile,
    float mins, float maxs ,float length, float thickness, float* offset, float t, int top) {

    int n_s_segments = 50;
    float s=mins;
    float delta_s = (maxs-mins)/((float) n_s_segments );

    //for (int i=0;i<n_s_segments;i++) {

        float x[24];  // a cube of coordinates
        float z[24];  // a cube of coordinates

        float nv[6];

//        float t =0.0;

        for(int n=0; n<n_s_segments;n++) {
            printf("floor flange segment %d %f - %f;  %f %f\n",n,mins,maxs, mins + (n)*delta_s,t);
            //for(int i=0;i<=1;i++) {
                int i=0;
                for(int j=0;j<=1;j++) {
                    for(int k=0;k<=1;k++) {
                        x[  12*i + 6*j + 3*k +0 ] = t +  (top?1.0:-1.0)*((float)j)*thickness ;
                        x[  12*i + 6*j + 3*k +1 ] = mins + (n+k)*delta_s ;
                        x[  12*i + 6*j + 3*k +2 ] = 0.0 ; //(((float)j))*length;
                     printf("x[%d] %f %d %d %d %d %f\n",12*i + 6*j + 3*k +2, x[  12*i + 6*j + 3*k +2 ] ,i,j,k,n,delta_s);
                     printf("%f %f %f\n",x[  12*i + 6*j + 3*k +0 ] ,x[  12*i + 6*j + 3*k +1 ] ,x[  12*i + 6*j + 3*k +2 ] );
                    }
                }
            //}

            for(int j=0;j<=1;j++) {
                unfaceted_surface_norm(  t  , mins + (n+j)*delta_s , spec, &nv[j*3] );
                xz_project_unit( &nv[j*3] );
                scale(&nv[j*3], length );
            }


            for(int i=0;i<4;i++) {
                faceted_X(&z[3*i],&x[3*i], spec ) ;
//                trnsfrm(&z[3*i],&x[3*i] ) ;
                add(&z[3*i],&z[3*i],offset);
                add( &z[12+ 3*i], &z[3*i], &nv[ 3*(i%2) ] );
            }

            if (top) {
                stl_triangle( &z[ 0 ], &z[12], &z[3] ,stlfile );  
                stl_triangle( &z[ 3 ], &z[12], &z[15] ,stlfile );  

                stl_triangle( &z[ 6 ], &z[9], &z[18] ,stlfile );  
                stl_triangle( &z[ 9 ], &z[21], &z[18] ,stlfile );  

                stl_triangle( &z[ 12 ], &z[18], &z[15] ,stlfile );  
                stl_triangle( &z[ 18 ], &z[21], &z[15] ,stlfile );  

                stl_triangle( &z[ 0 ], &z[3], &z[6] ,stlfile );  
                stl_triangle( &z[ 6 ], &z[3], &z[9] ,stlfile );  

                if (n==0) {
                    stl_triangle( &z[ 0 ], &z[6], &z[12] ,stlfile );  
                    stl_triangle( &z[ 6 ], &z[18], &z[12] ,stlfile );  
                }

                if (n==n_s_segments-1) {
                    printf("endcap\n");
                    stl_triangle( &z[ 3 ], &z[15], &z[9] ,stlfile );  
                    stl_triangle( &z[ 9 ], &z[15], &z[21] ,stlfile );  

                }

            } else {
                stl_triangle( &z[ 0 ], &z[3], &z[12] ,stlfile );  
                stl_triangle( &z[ 12 ], &z[3], &z[15] ,stlfile );  

                stl_triangle( &z[ 6 ], &z[18], &z[9] ,stlfile );  
                stl_triangle( &z[ 9 ], &z[18], &z[21] ,stlfile );  

                stl_triangle( &z[ 12 ], &z[15], &z[18] ,stlfile );  
                stl_triangle( &z[ 18 ], &z[15], &z[21] ,stlfile );  

                stl_triangle( &z[ 0 ], &z[6], &z[3] ,stlfile );  
                stl_triangle( &z[ 6 ], &z[9], &z[3] ,stlfile );  

                if (n==0) {
                    stl_triangle( &z[ 0 ], &z[12], &z[6] ,stlfile );  
                    stl_triangle( &z[ 6 ], &z[12], &z[18] ,stlfile );  
                }

                if (n==n_s_segments-1) {
                    printf("endcap\n");
                    stl_triangle( &z[ 3 ], &z[9], &z[15] ,stlfile );  
                    stl_triangle( &z[ 9 ], &z[21], &z[15] ,stlfile );  

                }


            }
/*
            printf("triangle",z[0],z[1],z[2]);
            stl_printquad( &z[6] , &z[0], stlfile, 0 ) ;
            stl_printquad( &z[12] , &z[18], stlfile, 0 ) ;

            stl_triangle( &z[ 0 ], &z[3], &z[12] ,stlfile );  //edge of the thing
            stl_triangle( &z[ 12 ], &z[3], &z[15] ,stlfile );  //edge of the thing

            stl_triangle( &z[ 6 ], &z[18], &z[9] ,stlfile );  //edge of the thing
            stl_triangle( &z[ 18 ], &z[21], &z[9] ,stlfile );  //edge of the thing
            */

        }


}


void write_texture_back_stl( void(*trnsfrm)(float*, float*), float mins, float maxs , float thickness, float* offset, stl_output_config cf) {

    int n_s_segments = 50;
    int n_t_segments = 50;

    float delta_s = (maxs-mins)/((float) n_s_segments );
    float delta_t = 1.0/((float) n_t_segments );

    float s=mins;
    float t=0.0;

    float x[24];  // a cube of coordinates
    float z[24];  // a cube of coordinates

    for(int n=0; n<n_s_segments;n++) {
    for(int m=0; m<n_t_segments;m++) {
      //  printf( "nm %d %d\n", n,m );
    
        for(int i=0;i<=1;i++) {
            for(int j=0;j<=1;j++) {
                for(int k=0;k<=1;k++) {
                    x[  12*i + 6*j + 3*k +0 ] = (m+k)*delta_t ;
                    x[  12*i + 6*j + 3*k +1 ] = mins + (n+j)*delta_s ;
                    x[  12*i + 6*j + 3*k +2 ] = -(((float)i))*thickness;
                }
            }
        }
    

        for(int i=0;i<8;i++) {
        	trnsfrm(&z[3*i],&x[3*i] ) ;
            add(&z[3*i],&z[3*i],offset);
        }

       // stl_printquad( &z[0] , &z[6], stlfile, 0 ) ;
        output_quad( &z[18] , &z[12], cf, 0 ) ;
        //stl_printquad( &x[12] , &x[18], stlfile, 1 ) ;

    } 
    }

}



void write_surface_stl( void(*trnsfrm)(float*, float*), FILE* stlfile,float mins, float maxs , float thickness, float* offset) {

    int n_s_segments = 50;
    int n_t_segments = 50;

    float delta_s = (maxs-mins)/((float) n_s_segments );
    float delta_t = 1.0/((float) n_t_segments );

    float s=mins;
    float t=0.0;

    float x[24];  // a cube of coordinates
    float z[24];  // a cube of coordinates

    for(int n=0; n<n_s_segments;n++) {
    for(int m=0; m<n_t_segments;m++) {
        //printf( "nm %d %d\n", n,m );
        
            for(int i=0;i<=1;i++) {
                for(int j=0;j<=1;j++) {
                    for(int k=0;k<=1;k++) {
                        x[  12*i + 6*j + 3*k +0 ] = (m+k)*delta_t ;
                        x[  12*i + 6*j + 3*k +1 ] = mins + (n+j)*delta_s ;
                        x[  12*i + 6*j + 3*k +2 ] = -(((float)i))*thickness;
                      //  printf("x[%d] %f %d %d %d %d %d %f %f\n",12*i + 6*j + 3*k +2, x[  12*i + 6*j + 3*k +2 ] ,i,j,k,n,m,delta_t,delta_s);
                      //  printf("%f %f %f\n",x[  12*i + 6*j + 3*k +0 ] ,x[  12*i + 6*j + 3*k +1 ] ,x[  12*i + 6*j + 3*k +2 ] );
                    }
                }
            }
        

        for(int i=0;i<8;i++) {
        	trnsfrm(&z[3*i],&x[3*i] ) ;
            add(&z[3*i],&z[3*i],offset);
        }

        stl_printquad( &z[0] , &z[6], stlfile, 0 ) ;
        stl_printquad( &z[18] , &z[12], stlfile, 0 ) ;
        //stl_printquad( &x[12] , &x[18], stlfile, 1 ) ;

        if (n==0) {
            stl_triangle( &z[ 0 ], &z[12], &z[3] ,stlfile );  //edge of the thing
            stl_triangle( &z[ 12 ], &z[15], &z[3] ,stlfile );  //edge of the thing
            //stl_triangle( &ray2[3], &ray2_back[3], &ray1_back[3],stlfile );  //edge of the thing
        }
        if (n==n_s_segments-1) {
            stl_triangle( &z[ 6 ], &z[9], &z[18] ,stlfile );  //edge of the thing
            stl_triangle( &z[ 18 ], &z[9], &z[21] ,stlfile );  //edge of the thing
            //stl_triangle( &ray2[3], &ray2_back[3], &ray1_back[3],stlfile );  //edge of the thing
        }

    } 
    }

}


void mold_parting_surface_norm(float u, int quadrant, pottoy_spec_t* spec, float* result) {


    quadrant = (quadrant + 4)%4;
    float s1=0; //= (quadrant & 1) ? -1.0 : 1.0 ;
    float s2=0; //= (quadrant & 2) ? 1.0 : -1.0 ;

    switch (quadrant) {
        case 0:
            s1 = 1.0;
            s2 = 1.0;

            break;
        case 1:
            s1 = -1.0;
            s2 = 1.0;

            break;
        case 2:
            s1 = -1.0;
            s2 = -1.0;

            break;
        case 3:
            s1 = 1.0;
            s2 = -1.0;
            break;
    }


    float curve_p[2];

//	spline_getPoint( uv[0], &point[0], spec );
    catmullrom2(u,&spec->points[0], spec->npoints, &curve_p[0]);
    
	float r = curve_p[0];
	float y = curve_p[1];

    float b_squared = spec->squash*spec->squash;
    float cos_theta = 1.0 / sqrt(1.0 + b_squared) ;
    float sin_theta = spec->squash / sqrt(1.0 + b_squared) ;

//    v1[0] = s1 * r * cos_theta ;
//    v1[1] = y;
//    v1[2] = s2 * r * spec->squash*sin_theta;

    //find the surface normal.
    float t1[3];  // 45 degrees to x and y
    float t2rs[2];  // along the spline
    float t2[3];  // along the spline

    t1[0] = -s2;
    t1[1] =  0.0;
    t1[2] =  s1;

    catmullrom2_tangent(u, &spec->points[0], spec->npoints, &t2rs[0]);

    //printf("spline t %f %f %f %f %f\n",u, curve_p[0], curve_p[1], t2rs[0], t2rs[1]);
    fflush(stdout);
    t2[0] = s1 * t2rs[0] * cos_theta ;
    t2[1] = t2rs[1];
    t2[2] = s2 * t2rs[0] * spec->squash*sin_theta;

    float nn[3];
    float n[3];
    float rr[3];

    cross_product(&t2[0],&t1[0],&nn[0]);
    normalize(&nn[0],&n[0]);

    cross_product(&t2[0],&n[0],&rr[0]);
    normalize(rr,result);

//    printf("norm %f %f %f %f  -- %f\n",u,n[0], n[1], n[2], norm(&n[0]) );
}


// this is the unfaceted surface only
void mold_parting_norm_ray( float u, float l,float thickness, int quadrant  , pottoy_spec_t* spec, float* result) {

    quadrant = (quadrant + 4)%4;
    float s1=0; //= (quadrant & 1) ? -1.0 : 1.0 ;
    float s2=0; //= (quadrant & 2) ? 1.0 : -1.0 ;

    switch (quadrant) {
        case 0:
            s1 = 1.0;
            s2 = 1.0;
//            s1 = 1.0;
 //           s2 = -1.0;
            break;
        case 1:
            s1 = -1.0;
            s2 = 1.0;

            break;
        case 2:
            s1 = -1.0;
            s2 = -1.0;

            break;
        case 3:
            s1 = 1.0;
            s2 = -1.0;
            break;
    }

    // on the surface, at points where the tangent plane makes a 45 degree angle with the x and y axes.
    float v1[3];

    // standing out along the normal from v1 and v2
    float v2[3];

    float curve_p[2];

//	spline_getPoint( uv[0], &point[0], spec );
    catmullrom2(u,&spec->points[0], spec->npoints, &curve_p[0]);
    
	float r = curve_p[0];
	float y = curve_p[1];

    float b_squared = spec->squash*spec->squash;
    float cos_theta = 1.0 / sqrt(1.0 + b_squared) ;
    float sin_theta = spec->squash / sqrt(1.0 + b_squared) ;

    v1[0] = s1 * r * cos_theta ;
    v1[1] = y;
    v1[2] = s2 * r * spec->squash*sin_theta;

    //find the surface normal.
    float t1[3];  // 45 degrees to x and y
    float t2rs[2];  // along the spline
    float t2[3];  // along the spline

    t1[0] = -s2;
    t1[1] =  0.0;
    t1[2] =  s1;

    catmullrom2_tangent(u, &spec->points[0], spec->npoints, &t2rs[0]);

    //printf("spline t %f %f %f %f %f\n",u, curve_p[0], curve_p[1], t2rs[0], t2rs[1]);
    fflush(stdout);
    t2[0] = s1 * t2rs[0] * cos_theta ;
    t2[1] = t2rs[1];
    t2[2] = s2 * t2rs[0] * spec->squash*sin_theta;

    float nn[3];
    float n[3];

    cross_product(&t2[0],&t1[0],&nn[0]);
    normalize(&nn[0],&n[0]);
    printf("norm %f %f %f %f  -- %f\n",u,n[0], n[1], n[2], norm(&n[0]) );

    v2[0] = v1[0] + l*n[0];
    v2[1] = v1[1] + l*n[1];
    v2[2] = v1[2] + l*n[2];

    v1[0] = v1[0] - thickness*n[0];
    v1[1] = v1[1] - thickness*n[1];
    v1[2] = v1[2] - thickness*n[2];



    result[0] = v1[0]; 
    result[1] = v1[1]; 
    result[2] = v1[2]; 
    result[3] = v2[0]; 
    result[4] = v2[1]; 
    result[5] = v2[2]; 

}

void faceted_X(float* xyz, float* uv, pottoy_spec_t* spec ) {

    if (spec->facet == 0) {
        sor( xyz, uv, spec );
        return;
    } 

    float point1[3];
    float point2[3];
    float point3[3];

    float uvA[2];
    float uvB[2];
    float uvC[2];

    float u = uv[0];
    float v0 = uv[1]; // +1.0; // the +1 here tries to fix phase problems with negative values of v
    float pufff =   (1.0-exp(-u*20.0))*(1-exp(-1.0*(1.0-u)*20.0)); //min(1,10*u,10*(1-u)) ;
	float v = v0 + u * (spec->twist) ;

    double ip;
	float phase = modf( ( (v0) * ((float) spec->n_facets)),&ip);
    if (phase<0) {phase += 1.0; } 
    float deltav = phase / ((float) spec->n_facets);
	float va = v-deltav;
    float vb = va+(1.0/((float) spec->n_facets));
	float f = phase ; // ((v - va)%1) * n_facets; 

    uvA[0] = u ; uvA[1] = va;
    uvB[0] = u ; uvB[1] = vb;
    uvC[0] = u ; uvC[1] = v ;

    sor( &point1[0], &uvA[0], spec );
    sor( &point2[0], &uvB[0], spec );
    sor( &point3[0], &uvC[0], spec );

    xyz[0] = point3[0] - pufff*(spec->puff-1.0)*( point1[0] + f*( point2[0] - point1[0] ) - point3[0]);
	xyz[1] = point3[1] - pufff*(spec->puff-1.0)*( point1[1] + f*( point2[1] - point1[1] ) - point3[1]);
	xyz[2] = point3[2] - pufff*(spec->puff-1.0)*( point1[2] + f*( point2[2] - point1[2] ) - point3[2]);

}


void surface_obj_bbox(float* bbox, pottoy_spec_t* spec) {
    float xyz[3];
    float uv[2];

    int n_sectors=90;
    int n_levels=50;
    int j=0;
    int i=0;
    double v = (((double) j)/((double) (n_sectors-1)));
    double u = ((double)i)/((double) n_levels);
    uv[0] = (float) u;
    uv[1] = (float) v;
    faceted_X(&xyz[0],&uv[0], spec);

    bbox[0] = xyz[0];
    bbox[1] = xyz[1];
    bbox[2] = xyz[2];

    bbox[3] = xyz[0];
    bbox[4] = xyz[1];
    bbox[5] = xyz[2];

    for(j=0;j<n_sectors;j++) {
        double v = (((double) j)/((double) (n_sectors-1)));
        //printf("v: %f\n",v);
        for(i=0;i<=n_levels;i++) {
            double u = ((double)i)/((double) n_levels);
            //double theta = u* 2.0*M_PI ;
            uv[0] = (float) u;
            uv[1] = (float) v;
            faceted_X(&xyz[0],&uv[0], spec);

            if (xyz[0] < bbox[0]) {bbox[0] = xyz[0];}
            if (xyz[1] < bbox[1]) {bbox[1] = xyz[1];}
            if (xyz[2] < bbox[2]) {bbox[2] = xyz[2];}
            if (xyz[0] > bbox[3]) {bbox[3] = xyz[0];}
            if (xyz[1] > bbox[4]) {bbox[4] = xyz[1];}
            if (xyz[2] > bbox[5]) {bbox[5] = xyz[2];}

        }
    }
}



// output an obj file augmented with some vertices having fixed UV coordiantes based on the parametric fn X()
int write_surface_obj2(char* fn, pottoy_spec_t* spec, int n_sectors, int n_levels) {

    FILE* objf = fopen(fn, "wt");
    float xyz[3];
    float uv[2];

//    int n_sectors=90;
//    int n_levels=50;

    for(int j=0;j<n_sectors;j++) {
        double v = (((double) j)/((double) (n_sectors-1)));
        //printf("v: %f\n",v);
        for(int i=0;i<=n_levels;i++) {
            double u = ((double)i)/((double) n_levels);
            //double theta = u* 2.0*M_PI ;
            uv[0] = (float) u;
            uv[1] = (float) v;
//            vt[1] = (float) theta;
            faceted_X(&xyz[0],&uv[0], spec);

        //void Xs(float* xyz, float* vtheta, pottoy_spec_t* spec );

//            Xs(&xyz[0], );

#if 1
            

// Pin the bottom corners to uv=(0,0) and (0,1)
// Pin the two sides of seam at theta=0.0 to u=0 and u=1 
            if (i==0 && j==0)  { // && (j==0 || j==n_levels-1)) {
                fprintf(objf, "v %f %f %f fix 0.0 0.0\n",xyz[0],xyz[1],xyz[2]);
            } else if (i==0 && j==n_sectors-1) {
                fprintf(objf, "v %f %f %f fix 1.0 0.0\n",xyz[0],xyz[1],xyz[2]);
            } else if (j==0 || j==n_sectors-1) {
                fprintf(objf, "v %f %f %f zip %d\n",xyz[0],xyz[1],xyz[2],i);
          //  } else if (i==n_sectors) { // && (j==0 || j==n_levels-1)) {
           //     fprintf(objf, "v %f %f %f fixu %f %d\n",xyz[0],xyz[1],xyz[2],u,j);  
            } else {
                fprintf(objf, "v %f %f %f\n",xyz[0],xyz[1],xyz[2]);
            }

            //stash the parameters v,theta in the place of texture coordinates
            fprintf(objf, "vn %f %f 0.0\n",u,v);
#else
        fprintf(objf,"v %f %f %f\n",xyz[0],xyz[1],xyz[2]);
#endif
        }
    }

    for(int j=0;j<n_sectors-1;j++) {
        for(int i=0;i<n_levels;i++) {
            fprintf(objf, "f %d %d %d\n",j*(n_levels+1)+i+1 , j*(n_levels+1)+i+2, (j+1)*(n_levels+1)+i+1);
            fprintf(objf, "f %d %d %d\n",i+1+(j+1)*(n_levels+1) , j*(n_levels+1)+i+2,i+2+(j+1)*(n_levels+1));
        }  
    }
    fclose(objf);
    return(0);
}

// output an obj file augmented with some vertices having fixed UV coordiantes based on the parametric fn X()
int write_surface_obj0(char* fn, pottoy_spec_t* spec, int n_sectors, int n_levels) {

    FILE* objf = fopen(fn, "wt");
    float xyz[3];
    float uv[2];

//    int n_sectors=90;
//    int n_levels=50;

    for(int j=0;j<n_sectors;j++) {
        double v = (((double) j)/((double) (n_sectors-1)));
        //printf("v: %f\n",v);
        for(int i=0;i<=n_levels;i++) {
            double u = ((double)i)/((double) n_levels);
            //double theta = u* 2.0*M_PI ;
            uv[0] = (float) u;
            uv[1] = (float) v;
            
//            vt[1] = (float) theta;
            faceted_X(&xyz[0],&uv[0], spec);
    //    printf("%f %f -> %f %f %f\n",uv[0],uv[1],xyz[0],xyz[1],xyz[2]);
        //void Xs(float* xyz, float* vtheta, pottoy_spec_t* spec );

//            Xs(&xyz[0], );

        fprintf(objf,"v %f %f %f\n",xyz[0],xyz[1],xyz[2]);
        }
    }

    for(int j=0;j<n_sectors-1;j++) {
        for(int i=0;i<n_levels;i++) {
            fprintf(objf, "f %d %d %d\n",j*(n_levels+1)+i+1 , j*(n_levels+1)+i+2, (j+1)*(n_levels+1)+i+1);
            fprintf(objf, "f %d %d %d\n",i+1+(j+1)*(n_levels+1) , j*(n_levels+1)+i+2,i+2+(j+1)*(n_levels+1));
        }  
    }
    fclose(objf);
    return(0);
}





// output an obj file augmented with some vertices having fixed UV coordiantes based on the parametric fn X()
int write_surface_obj(char* fn) {

    FILE* objf = fopen(fn, "wt");
    float xyz[3];
    float vt[2];

    int n_sectors=90;
    int n_levels=50;

    for(int j=0;j<n_levels;j++) {
        double v = (((double) j)/((double) (n_levels-1)));

        for(int i=0;i<=n_sectors;i++) {
            double u = ((double)i)/((double) n_sectors);
            double theta = u* 2.0*M_PI ;
            vt[0] = (float) v;
            vt[1] = (float) theta;
            X(&xyz[0],&vt[0]);

#if 1
            

// Pin the bottom corners to uv=(0,0) and (0,1)
// Pin the two sides of seam at theta=0.0 to u=0 and u=1 
            if (i==0 && j==0)  { // && (j==0 || j==n_levels-1)) {
                fprintf(objf, "v %f %f %f fix 0.0 0.0\n",xyz[0],xyz[1],xyz[2]);
            } else if (i==n_sectors && j==0) {
                fprintf(objf, "v %f %f %f fix 1.0 0.0\n",xyz[0],xyz[1],xyz[2]);
            } else if (i==0) {
                fprintf(objf, "v %f %f %f fixu %f %d\n",xyz[0],xyz[1],xyz[2],u,j);
            } else if (i==n_sectors) { // && (j==0 || j==n_levels-1)) {
                fprintf(objf, "v %f %f %f fixu %f %d\n",xyz[0],xyz[1],xyz[2],u,j);  
            } else {
                fprintf(objf, "v %f %f %f\n",xyz[0],xyz[1],xyz[2]);
            }

            //record the parameters v,theta in the place of texture coordinates
            fprintf(objf, "vt %f %f 0.0\n",v,theta);
#else
        fprintf(objf,"v %f %f %f\n",xyz[0],xyz[1],xyz[2]);
#endif
        }
    }

    for(int j=0;j<n_levels-1;j++) {
        for(int i=0;i<n_sectors;i++) {
            fprintf(objf, "f %d %d %d\n",j*(n_sectors+1)+i+1 , j*(n_sectors+1)+i+2, (j+1)*(n_sectors+1)+i+1);
            fprintf(objf, "f %d %d %d\n",i+1+(j+1)*(n_sectors+1) , j*(n_sectors+1)+i+2,i+2+(j+1)*(n_sectors+1));
        }  
    }
    fclose(objf);
    return(0);
}





static void scan_array(const char *str, int len, void *user_data) {
    struct json_token t;
    int i;
  //  printf("Parsing array: %.*s\n", len, str);
    float* result = (float*) user_data;
    int ir=0;
    for (i = 0; json_scanf_array_elem(str, len, "", i, &t) > 0; i++) {

        float x;
        float y;
        json_scanf(t.ptr, t.len, "{ x:%f, y:%f }",
             &result[ir],  &result[ir+1]);
 //     printf("Index %d, token [%.*s] %f %f\n", i, t.len, t.ptr, result[ir], result[ir+1]);
        ir+=2;
    }
}

void scale_spec(float scale, pottoy_spec_t* spec) {
    for (int i = 0; i< spec->npoints*2; i++) {
        spec->points[i] *= scale;
    }

}

int read_spec(char* filename, pottoy_spec_t* spec ) {


    /* declare a file pointer */
    FILE    *infile;
    char    *buffer;
    long    numbytes;
    
    /* open an existing file for reading */
    infile = fopen( filename , "r");
    /* quit if the file does not exist */
    if(infile == NULL)
        return 1;

    //printf("opened file %s\n", filename);

    /* Get the number of bytes */
    fseek(infile, 0L, SEEK_END);
    numbytes = ftell(infile);
    printf("need %d bytes\n",numbytes);  fflush(stdout);
    /* reset the file position indicator to 
    the beginning of the file */
    fseek(infile, 0L, SEEK_SET);	
    
    /* grab sufficient memory for the 
    buffer to hold the text */
    buffer = (char*)calloc(numbytes+1, sizeof(char));	
    
    /* memory error */
    if(buffer == NULL)
        return 1;
    
    /* copy all the text into the buffer */
    int bytes_read = fread(buffer, sizeof(char), numbytes, infile);
    buffer[numbytes]=0;
    printf("read %d bytes\n",bytes_read);  fflush(stdout);
    fclose(infile);
    printf("buffer: %s\n",buffer); fflush(stdout);


  int a = 0, d = 0;
  int n_points = 0;
  char *c = NULL;
  void *my_data = NULL;
  float* points;

  //static const char str[] = "{ \"a\": 123, \"d\": false, \"b\": [1, 2], \"c\": \"hi\" } ";

  //  char* str = ;
  // printf("input: %s\n",buffer);
 
    json_scanf(buffer, strlen(buffer), "{ n_facets:%d, facet:%B, npoints:%d, puff:%f, squash:%f, twist:%f }",
             &spec->n_facets, &spec->facet, &spec->npoints,
            &spec->puff, &spec->squash, &spec->twist
             );

    printf("n_facets: %d\n",spec->n_facets);
    printf("facet: %d\n",spec->facet);
   printf("npoints: %d\n",spec->npoints);

    spec->points = (float*)calloc( spec->npoints, 2*sizeof(float) );


      json_scanf(buffer, strlen(buffer), "{ points:%M  }",
             scan_array, spec->points );

    for (int i = 0; i< n_points*2; i++) {
        printf("xy %f\n",points[i]);
    }

    free(buffer);
    return(0);

}

void add_triangle_contours(TESStesselator* tess, pottoy_spec_t* spec, int mn_sectors, int mn_levels) {

    float x[6];
    float uv[2];

#define UU(i) (((double)(i))/((double) mn_levels))
#define VV(j) (((double) (j))/((double) (mn_sectors-1)))

    for(int j=-mn_sectors;j<2*mn_sectors-1;j++) {
        for(int i=0;i<mn_levels;i++) {

            x[0] = UU(i);   x[1] = VV(j);
            x[2] = UU(i+1); x[3] = VV(j);
            x[4] = UU(i);   x[5] = VV(j+1);

            //printf("z %f,%f %f,%f %f,%f \n",x[0],x[1],x[2],x[3],x[4],x[5]);
            tessAddContour(tess, 2, (void*) &(x[0]), sizeof(float)*2, 3);

            x[0] = UU(i)  ; x[1] = VV(j+1);
            x[2] = UU(i+1); x[3] = VV(j);
            x[4] = UU(i+1); x[5] = VV(j+1);

            //printf("z %f,%f %f,%f %f,%f \n",x[0],x[1],x[2],x[3],x[4],x[5]);
            tessAddContour(tess, 2, (void*) &(x[0]), sizeof(float)*2, 3);

//            fprintf(objf, "f %d %d %d\n",j*(n_levels+1)+i+1 , j*(n_levels+1)+i+2, (j+1)*(n_levels+1)+i+1);
//            fprintf(objf, "f %d %d %d\n",i+1+(j+1)*(n_levels+1) , j*(n_levels+1)+i+2,i+2+(j+1)*(n_levels+1));
        }  
    }
    fflush(stdout);
}

 

int contour_to_mesa(int n, float* p, float* offset, float f, int nt, float* result_buffer) {
    float center[3] = {0.0,0.0,0.0};
    float c[3] = {0.0,0.0,0.0};
    float d[3] = {0.0,0.0,0.0};
    float o[3] = {0.0,0.0,0.0};

    float* a;
    float* b;
    assert(nt>=3*n);
//    memcpy(result_buffer,p,n*3*sizeof(float));

    for (int i=0;i<n;i++) {
        center[0] += p[i*3 + 0];
        center[1] += p[i*3 + 1];
        center[2] += p[i*3 + 2];
    }
    center[0] /= n;
    center[1] /= n; 
    center[2] /= n; 
   // set(&o[0],&center[0]);
    add(&o[0],&center[0],offset);
/*
 
      |a      c
    i *------*-------* o   center
      |        \    /
      |       ___>* 
      |  ___/      d
      |b/
  i+1 *
      |
    
    triangles   abc ,  cbd,  cdo

*/

    int nr = 0;
    for (int i=0;i<n;i++) {
        a = &p[i*3]        ;
        b = &p[((i+1)%n)*3];

        weighted_sum3( &c[0] , f, 1-f, a, &center[0]);	
        weighted_sum3( &d[0] , f, 1-f, b, &center[0]);	

        add( &c[0], &c[0], offset );
        add( &d[0], &d[0], offset );

        set_triangle( &result_buffer[nr*9], a,b,c);  nr++;
        set_triangle( &result_buffer[nr*9], c,b,d);  nr++;
        set_triangle( &result_buffer[nr*9], c,d,o);  nr++;

    }
    return(nr);

}
