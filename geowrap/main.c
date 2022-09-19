#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
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


void print_triangle_bbox( float* z1, float* z2, float* z3 , FILE* fp, float* bbox) {

	float v1[3];
	float v2[3];
	float v3[3];

	transform(v1,z1);
	transform(v2,z2);
	transform(v3,z3);

    float miny = bbox[1];
    float maxy = bbox[4];

    if (v1[1] > maxy && v2[1] > maxy && v3[1] > maxy  ) {return;}
    if (v1[1] < miny && v2[1] < miny  && v3[1] < miny ) {return;}

	float n[3];
	triangle_normal(v1,v2,v3,&n[0]);
//    printf("print triangle with normal %f %f %f\n",n[0], n[1], n[2] );
	if (n[0]==0.0 && n[1]==0.0 && n[2]==0.0) { 
        //printf("zero normla\n") ;; 
        return;}

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
	fprintf(fp,"      vertex %f %f %f\n",v1[0],v1[1],v1[2]);
	fprintf(fp,"      vertex %f %f %f\n",v2[0],v2[1],v2[2]);
	fprintf(fp,"      vertex %f %f %f\n",v3[0],v3[1],v3[2]);
	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");

}


int main(int argc, char *argv[])
{

	int n_sectors = 90;
	int n_levels = 50;

    if (!( argc==3 || argc==6)) {
        printf("usage: geowrap <input.json> <output.obj>\n");
        printf("usage: geowrap <texture_geometry.obj> <input-lscm.obj> <out.stl> <surf.json> <width>\n"); //<input.json> <output.obj>\n");
        exit(0);
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

    if (argc==6) {

        // argv[1] -- .obj file to be wrapped
        // argv[2] -- obj file with the uv coordinates from lscm
        // argv[3] -- output stl file
        // argv[4] -- surface spec in .json
        // argv[5] -- wrapping width of the texture
 
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
        MeshTriangles* mt = parse_triangles(argv[2],(float) width);
        printf("parsed interpolation\n");


        MeshTriangles* texture = parse_triangles_with_normals(argv[1],(float) 1.0);
        printf("parsed with normals\n");
        apply_interpolation_to_mesh( texture, mt, 10.0 );
        //apply_transform_to_mesh( texture, transf_X );
        //write_to_stl(texture,stlfile);


        for (int i=0; i<texture->ntriangles; i++) {
            //printf("mt path %d\n",i);
            int n1 = texture->triangles[3*i];
            int n2 = texture->triangles[3*i+1];
            int n3 = texture->triangles[3*i+2];
            //tessAddContour(tess, 2, &path_points[ path_offsets[i] ] , sizeof(float)*2, path_lengths[i]);
            float v1[3] = { texture->xpoints[n1*3+0] ,texture->xpoints[n1*3+1] ,texture->xpoints[n1*3+2] };
            float v3[3] = { texture->xpoints[n2*3+0] ,texture->xpoints[n2*3+1] ,texture->xpoints[n2*3+2]};
            float v2[3] = { texture->xpoints[n3*3+0] ,texture->xpoints[n3*3+1] ,texture->xpoints[n3*3+2] };

//            float v2[3] = { texture->paths[i*6+2] ,texture->paths[i*6+3] , 0.1 };
//            float v3[3] = { texture->paths[i*6+4] ,texture->paths[i*6+5] , 0.1 };
            print_triangle_bbox(v1,v2,v3, stlfile, bbox);
        }





        // print out the shell of the surface itself -- outside
        for (int i=0; i<mt->npaths; i++) {
            //printf("mt path %d\n",i);
            //tessAddContour(tess, 2, &path_points[ path_offsets[i] ] , sizeof(float)*2, path_lengths[i]);
            float v1[3] = { mt->paths[i*6+0] ,mt->paths[i*6+1] , -0.1 };
            float v2[3] = { mt->paths[i*6+2] ,mt->paths[i*6+3] , -0.1 };
            float v3[3] = { mt->paths[i*6+4] ,mt->paths[i*6+5] , -0.1 };
            print_triangle(v1,v2,v3, stlfile);
        }


	fprintf(stlfile, "endsolid x\n");
    fclose(stlfile);
    }

}

