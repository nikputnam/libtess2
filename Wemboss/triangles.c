
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "triangles.h"
#include "vectorops.h"
#include "meshindex2d.h"

void add_triangle(float* v1, float* v2, float* v3 , MeshTriangles* mt) {
    printf("add triangle %d %d / %d %d\n",mt->ntriangles,mt->nxpoints, mt->max_ntriangles, mt->max_nxpoints);
	while ((mt->nxpoints+3)>=mt->max_nxpoints) {
		mt->max_nxpoints *= 2;		printf("more points %d\n",mt->max_nxpoints);  fflush(stdout);
		mt->xpoints = (float*) realloc( mt->xpoints, sizeof(float)*mt->max_nxpoints*3);
	}

	while ((mt->ntriangles+1)>=mt->max_ntriangles) {
		mt->max_ntriangles *= 2;
		printf("more triangles %d\n",mt->max_ntriangles);
        fflush(stdout);
		mt->triangles = (float*) realloc( mt->triangles, sizeof(int)*mt->max_ntriangles*3);
	}

	mt->triangles[ (mt->ntriangles * 3)  ] = mt->nxpoints;
	mt->triangles[ (mt->ntriangles * 3)+1] = mt->nxpoints+1;
	mt->triangles[ (mt->ntriangles * 3)+2] = mt->nxpoints+2;
	mt->ntriangles++;

//    return;

	set( &mt->xpoints[ (mt->nxpoints) * 3  ], v1 );
 	set( &mt->xpoints[ (mt->nxpoints+1) * 3  ], v2 );
	set( &mt->xpoints[ (mt->nxpoints+2) * 3  ], v3 );

	mt->nxpoints+=3;
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
	if (l==0.0) { 
        printf("zero area triangle?\n");
		n[0]=0.0;
		n[1]=0.0;
		n[2]=0.0;
		return;
	}
	n[0] /= l;
	n[1] /= l;
	n[2] /= l;
}


void print_triangle_raw( float* v1, float* v2, float* v3 , FILE* fp) {


	float n[3];
	triangle_normal(v1,v2,v3,&n[0]);

	//printf("  facet normal %f %f %f\n", n[0], n[1], n[2]);
	//fprintf(fp,"  facet normal %f %f %f\n", 0.0,0.0,0.0);
    fprintf(fp,"  facet normal %f %f %f\n", n[0], n[1], n[2]);
	fprintf(fp,"    outer loop\n");
	fprintf(fp,"      vertex %f %f %f\n",v1[0],v1[1],v1[2]);
	fprintf(fp,"      vertex %f %f %f\n",v2[0],v2[1],v2[2]);
	fprintf(fp,"      vertex %f %f %f\n",v3[0],v3[1],v3[2]);
	fprintf(fp,"    endloop\n");
	fprintf(fp,"  endfacet\n");


}

void apply_transform_to_mesh(MeshTriangles*  t,  void(*trnsfrm)(float*, float*) ) {
    for (int i=0;i<t->nxpoints;i++) {

        printf("%f %f %f\n", t->xpoints[3*i], t->xpoints[3*i+1], t->xpoints[3*i+2]);
        //t->xpoints[3*i+2] = 0.0;
        trnsfrm( & t->xpoints[3*i], & t->xpoints[3*i] );
        printf("                                   -> %f %f %f\n", t->xpoints[3*i], t->xpoints[3*i+1], t->xpoints[3*i+2]);

        if (i%10000 == 0) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d / %d = %f",i,t->nxpoints, ((float)i) / t->nxpoints) ;   fflush(stdout);}
        printf("\n");
    }

}

void apply_interpolation_to_mesh( MeshTriangles* t, MeshTriangles* interpolater, float zscale, meshindex* mi ) {

#if 1

    char* flags;
    flags = calloc( t->nxpoints, sizeof(char) );

    for (int i=0;i<t->nxpoints;i++) {
        flags[1]=0;
    }

    int n_hit = 0;
    int n_tries = 0;

    for (int i=0;i<t->ntriangles;i++) {
//        fprintf(stlfile,"triangle %d\n",i);
        int v1 = t->triangles[i*3];
        int v2 = t->triangles[i*3+1];
        int v3 = t->triangles[i*3+2];

        if (!flags[v1]) { flags[v1]=1;  mesh_interpolation(interpolater, & t->xpoints[3*v1], & t->xpoints[3*v1],mi, &n_hit ); n_tries++;}
        if (!flags[v2]) { flags[v2]=1;  mesh_interpolation(interpolater, & t->xpoints[3*v2], & t->xpoints[3*v2],mi, &n_hit ); n_tries++;}
        if (!flags[v3]) { flags[v3]=1;  mesh_interpolation(interpolater, & t->xpoints[3*v3], & t->xpoints[3*v3],mi, &n_hit ); n_tries++;}

    }

    printf("hit count: %d / %d\n",n_hit, n_tries);
#else
    for (int i=0;i<t->nxpoints;i++) {

       // printf("%f %f %f\n", t->xpoints[3*i], t->xpoints[3*i+1], t->xpoints[3*i+2]);
        mesh_interpolation(interpolater, & t->xpoints[3*i], & t->xpoints[3*i] );
        t->xpoints[3*i+2] = zscale*t->xpoints[3*i+2]; 
        //t->xpoints[3*i+1] = t->xpoints[3*i+1] * -1.0;
        //printf("                                   -> %f %f %f\n", t->xpoints[3*i], t->xpoints[3*i+1], t->xpoints[3*i+2]);

        if (i%10000 == 0) { printf("%d / %d = %f\n",i,t->nxpoints, ((float)i) / t->nxpoints) ; }
    }
#endif
}

void write_to_obj(MeshTriangles* t, FILE* stlfile, int yup) {

    if (yup) {

        for (int i=0;i<t->nxpoints;i++) {
            fprintf(stlfile,"v %f %f %f\n",t->xpoints[i*3 + 0],t->xpoints[i*3 + 2],-t->xpoints[i*3 + 1]);   
        }

        for (int i=0;i<t->ntriangles;i++) {
            fprintf(stlfile,"f %d %d %d\n",1+t->triangles[i*3 + 0],1+t->triangles[i*3 + 1],1+t->triangles[i*3 + 2]);   
        }

    } else {
        for (int i=0;i<t->nxpoints;i++) {
            fprintf(stlfile,"v %f %f %f\n",t->xpoints[i*3 + 0],t->xpoints[i*3 + 1],t->xpoints[i*3 + 2]);   
        }

        for (int i=0;i<t->ntriangles;i++) {
            fprintf(stlfile,"f %d %d %d\n",1+t->triangles[i*3 + 0],1+t->triangles[i*3 + 1],1+t->triangles[i*3 + 2]);   
        }
    }
}

void write_to_stl( MeshTriangles* t, FILE* stlfile ) {
    //FILE* stlfile = fopen(filename, "wt"); 
	//fprintf(stlfile,"solid x\n");

    for (int i=0;i<t->ntriangles;i++) {
//        fprintf(stlfile,"triangle %d\n",i);
        int v1 = t->triangles[i*3];
        int v2 = t->triangles[i*3+1];
        int v3 = t->triangles[i*3+2];
/*
        fprintf(stlfile,"%d %d %d\n",v1,v2,v3 );
        fprintf(stlfile,"%f %f %f\n",t->xpoints[v1*3], t->xpoints[v1*3+1],  t->xpoints[v1*3+2] );
        fprintf(stlfile,"%f %f %f\n",t->xpoints[v2*3], t->xpoints[v2*3+1],  t->xpoints[v2*3+2] );
        fprintf(stlfile,"%f %f %f\n",t->xpoints[v3*3], t->xpoints[v3*3+1],  t->xpoints[v3*3+2] );
*/
        print_triangle_raw( &t->xpoints[v1*3], &t->xpoints[v2*3], &t->xpoints[v3*3], stlfile);

    }

//	fprintf(stlfile, "endsolid x\n");
//    fclose(stlfile);
}

MeshTriangles* parse_triangles(char* filename, float width) {
    return parse_triangles_internal( filename, width, 0, 1 );
}

MeshTriangles* parse_triangles_with_normals(char* filename, float width) {
    return parse_triangles_internal( filename, width, 1, 0 );
}

void free_meshtriangles( MeshTriangles* mt ) {
    free(mt->points);
    free(mt->tpoints);
    free(mt->xpoints);
    free(mt->triangles);
    free(mt->ttriangles);
    free(mt->paths); 
    free(mt);
}

float mesh_width(MeshTriangles* mesh) {
//    float w=0;
    float minx = mesh->xpoints[0];
    float maxx = mesh->xpoints[0];


    for (int i=0; i<mesh->ntriangles; i++) {
        int v1 = mesh->triangles[i*3];
        int v2 = mesh->triangles[i*3+1];
        int v3 = mesh->triangles[i*3+2];

        float* x1 = &mesh->xpoints[v1*3];
        float* x2 = &mesh->xpoints[v2*3];
        float* x3 = &mesh->xpoints[v3*3];

        if (minx > x1[0]) {minx = x1[0];} 
        if (minx > x2[0]) {minx = x2[0];} 
        if (minx > x3[0]) {minx = x3[0];} 
        if (maxx < x1[0]) {maxx = x1[0];} 
        if (maxx < x2[0]) {maxx = x2[0];} 
        if (maxx < x3[0]) {maxx = x3[0];} 

    }

    return (maxx-minx);
}


MeshTriangles* parse_triangles_internal(char* filename, float width, int with_normals, int reduplicate) {

    MeshTriangles* mt = malloc(sizeof(MeshTriangles));

    FILE* file = fopen(filename, "r"); 
    if (file == NULL) {
        perror("Failed: ");
        printf("file: %s\n",filename); 
        //return 1;
        return (MeshTriangles* )NULL;
    }

    char line[1024];

    printf("file: %s\n",filename); fflush(stdout);
    //printf("%x\n",file);

    int npoints = 0;
    int ntriangles = 0;
    int ntpoints = 0;
    int nxpoints = 0;
    while (fgets(line, sizeof(line), file)) {
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
        if (strncmp(line,"vn ",3)==0) { 
            //printf("%s", line); 
            ntpoints++;
        }

        if (strncmp(line,"v ",2)==0) { 
            //printf("%s", line); 
            nxpoints++;
        }

        if (strncmp(line,"vt ",3)==0) { 
            //printf("%s", line); 
            npoints++;
        }

        if (strncmp(line,"f ",2)==0) { 
            //printf("%s", line);  fflush(stdout);
            ntriangles++;
        }
    }

    if (reduplicate) {
    // x3 because we're replicating the interpolation mesh 
        mt->npoints = 3*npoints;
        mt->ntpoints = 3*ntpoints;
        mt->nxpoints = nxpoints;
        mt->ntriangles = 3*ntriangles;
        mt->npaths = 3*ntriangles;
    } else {
        mt->npoints = npoints;
        mt->ntpoints = ntpoints;
        mt->nxpoints = nxpoints;
        mt->ntriangles = ntriangles;
        mt->npaths = ntriangles;

    }

    printf("npoints: %d\n",npoints);
    printf("ntpoints: %d\n",ntpoints);
    printf("nxpoints: %d\n",nxpoints);
    printf("ntriangles: %d\n",ntriangles);
    fflush(stdout);

    fclose(file);

    mt->points = malloc( sizeof(float)*(mt->npoints)*2  );  assert(mt->points != 0);
    mt->tpoints = malloc( sizeof(float)*(mt->ntpoints)*2  );assert(mt->tpoints != 0);
    mt->xpoints = malloc( sizeof(float)*(mt->nxpoints)*3 );assert(mt->xpoints != 0);
    mt->triangles = malloc( sizeof(int)*(mt->ntriangles)*3  );assert(mt->triangles != 0);
    mt->ttriangles = malloc( sizeof(int)*(mt->ntriangles)*3  );assert(mt->ttriangles != 0);
  //  mt->xtriangles = malloc( sizeof(int)*ntriangles*3 );assert(mt->xtriangles != 0);
    mt->paths = malloc( sizeof(float)*(mt->ntriangles)*6 );assert(mt->paths != 0);

    file = fopen(filename, "r");   assert(file);
    printf("x\n"); fflush(stdout);
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    int m=0;
    int xi=0;
    int maxtp = 0;
    while (fgets(line, sizeof(line), file)) {
                  //  printf("%d %d %d %d %s", i, l, j , m, line);  fflush(stdout);

        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
        if (strncmp(line,"vt",2)==0) { 
            float x;
            float y;
            //printf("%s", line); 
            sscanf(line,"vt %f %f",&x,&y);
            //printf("%f %f %s\n",x,y,line);
            mt->points[i++]=x*width;
            mt->points[i++]=-1.0*y*width;
//            mt->points[i++]=x;
//            mt->points[i++]=y;
        }

        if (strncmp(line,"v ",2)==0) { 
            float x;
            float y;
            float z;
            //printf("%s", line); 
            sscanf(line,"v %f %f %f",&x,&y,&z);
           // printf("l= %d ntpoints= %d %f %f %s\n",l,ntpoints,x,y,line);fflush(stdout);
            assert(xi<3*nxpoints);
            mt->xpoints[xi++]=x;
            mt->xpoints[xi++]=-1.0*z;
            mt->xpoints[xi++]=y;
            //printf("vv1: %f %f %f\n",x,-1.0*z,y);
        }


        if (strncmp(line,"vn ",3)==0) { 
            float x;
            float y;
            //printf("%s", line); 
            sscanf(line,"vn %f %f",&x,&y);
           // printf("l= %d ntpoints= %d %f %f %s\n",l,ntpoints,x,y,line);fflush(stdout);
            assert(l<2*ntpoints);
            mt->tpoints[l++]=x;
            mt->tpoints[l++]=y;
        }


        if (strncmp(line,"f ",2)==0) { 
            int a,b,c,x,y,z, n1,n2,n3;
            
            if ( 9==sscanf(line,"f %d/%d/%d %d/%d/%d %d/%d/%d",&a,&x,&n1,&b,&y,&n2,&c,&z,&n3) ) {
                
            mt->ttriangles[m++]=x-1;
            mt->ttriangles[m++]=y-1;
            mt->ttriangles[m++]=z-1;

            } else if (  
                6==sscanf(line,"f %d/%d %d/%d %d/%d",&a,&x,&b,&y,&c,&z)
            ) {

            mt->ttriangles[m++]=x-1;
            mt->ttriangles[m++]=y-1;
            mt->ttriangles[m++]=z-1;

            } else if (
                3==sscanf(line,"f %d %d %d",&a,&b,&c)

            ) {

            }

            mt->triangles[j++]=a-1;
            mt->triangles[j++]=b-1;
            mt->triangles[j++]=c-1;

            maxtp = maxtp < a-1 ? a-1 : maxtp;
            maxtp = maxtp < b-1 ? b-1 : maxtp;
            maxtp = maxtp < c-1 ? c-1 : maxtp;



        }
    }


    if (reduplicate) {
   
        //replicate image points two more times
        int ni = i;
        for ( int jj=0; jj<ni; jj+=2 ) { 
            mt->points[i++]=mt->points[jj  ] - width ;// +1.0; 
            mt->points[i++]=mt->points[jj+1];// - width;  

          //  printf("ttm: %f %f %d %d\n",mt->points[jj  ] - width, mt->points[jj+1], i,jj);

        }

        for ( int jj=0; jj<ni; jj+=2 ) { 
            mt->points[i++]=mt->points[jj  ] + width; // -1.0; 
            mt->points[i++]=mt->points[jj+1] ;//+ width;  

           // printf("ttm: %f %f %d %d\n",mt->points[jj  ] + width, mt->points[jj+1],i,jj);

        }
        assert(i==npoints*2 *3);

        //replicate image space triangles two more times
        int nj = j;
        for (int jj=0; jj<nj; jj++) {
            mt->triangles[j++]=mt->triangles[jj]+ni/2;
        } 

        for (int jj=0; jj<nj; jj++) {
            mt->triangles[j++]=mt->triangles[jj]+ni;
        } 
        assert(j==ntriangles*3 *3);

    

/*

    printf("ntriangles: %d %d %d %f\n",mt->ntriangles, nj,j, width);
    for (int i=0; i<mt->ntriangles;i++ ) {
        float* a = & mt->points[ mt->triangles[i*3]*2    ];
        float* b = & mt->points[ mt->triangles[i*3+1]*2 ];
        float* c = & mt->points[ mt->triangles[i*3+2]*2 ];
        printf("pid %d %d %d %d %d\n",i,0,mt->triangles[i*3]*2,npoints*2 *3, ni);
        printf("ttt %f %f\nttt %f %f\nttt %f %f\nttt\nttt\nttt\n",a[0],a[1], b[0],b[1], c[0],c[1]);

    }

//    exit(0);
*/
        //replicate texture points two more times
        int nl = l;
        for ( int jj=0; jj<nl; jj+=2 ) { 
            mt->tpoints[l++]=mt->tpoints[jj  ];//+1.0; 
            mt->tpoints[l++]=mt->tpoints[jj+1]-1.0;  
        }

        for ( int jj=0; jj<nl; jj+=2 ) { 
            mt->tpoints[l++]=mt->tpoints[jj  ];//-1.0; 
            mt->tpoints[l++]=mt->tpoints[jj+1]+1.0;  
        }


        //replicate texture space triangles two more times
        int nm = m;
        for (int jj=0; jj<nm; jj++) {
            mt->ttriangles[m++]=mt->ttriangles[jj]+nl/2;
        } 

        
        for (int jj=0; jj<nm; jj++) {
            mt->ttriangles[m++]=mt->ttriangles[jj]+nl;
        } 

    }

    printf("j got to %d; max triangle point %d ; ntriangles %d\n",j,maxtp, mt->ntriangles); fflush(stdout);

    if (reduplicate) {

    k=0;
    for (int i=0; i<mt->ntriangles;i++ ) {
            mt->paths[k++] = mt->tpoints[ mt->triangles[i*3]*2  ];
            mt->paths[k++] = mt->tpoints[ mt->triangles[i*3]*2+1];
            mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+1]*2 ];
            mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+1]*2 +1];
            mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+2]*2 ];
            mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+2]*2 +1];
    }
    }
       // assert(k==ntriangles*6*3);

    printf("%d %d\n",i,j);
    fflush(stdout);
    fclose(file);



/*
    for (int i=0; i<mt->ntriangles;i++ ) {
        float* a = & mt->points[ mt->triangles[i*3]*2    ];
        float* b = & mt->points[ mt->triangles[i*3+1]*2 ];
        float* c = & mt->points[ mt->triangles[i*3+2]*2 ];

        printf("ttt %f %f; %f %f; %f %f; \n",a[0],a[1], b[0],b[1], c[0],c[1]);
    }   
*/
    return mt;

}



float min(float a, float b, float c  ) {
    if (a<=b && a<=c) {return a;} 
    if (b<=a && b<=c) {return b;} 
    return c;
}

float max(float a, float b, float c  ) {
    if (a>=b && a>=c) {return a;} 
    if (b>=a && b>=c) {return b;} 
    return c;
}

//#define DEBUG_MI 

void mesh_interpolation(MeshTriangles* mt, float* p, float* uv, meshindex* mi, int* hit_count) {

//    uv[0]=0.0;
//    uv[1]=0.0;
//    return;
 
    float best_mm = 1000000000.0;
    float best_u = 10.0;
    float best_v = 10.0;
    float best_w = 10.0;
    int best_triangle = -1;

#ifdef DEBUG_MI
    float savep[2] = {p[0], p[1]};
#endif

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
            
        }*/

  //  for (int i=0; i<mt->ntriangles;i++ ) {


        meshindex_it* it = index_iterator( p[0],p[1],mi );
        int i = it->t;
        while (!it->done) { 
            i = it->t;

        float* a = & mt->points[ mt->triangles[i*3]*2    ];
        float* b = & mt->points[ mt->triangles[i*3+1]*2 ];
        float* c = & mt->points[ mt->triangles[i*3+2]*2 ];


        float* x = & mt->tpoints[ mt->ttriangles[i*3]*2    ];
        float* y = & mt->tpoints[ mt->ttriangles[i*3+1]*2 ];
        float* z = & mt->tpoints[ mt->ttriangles[i*3+2]*2 ];

        float v0[2];
        float v1[2];
        float v2[2];

        diff( c,a,&v0[0] );        
        diff( b,a,&v1[0] );        
        diff( p,a,&v2[0] );      

/*
        if ( (dist2(a,c) > 10.0) ||(dist2(b,a) > 10.0) || (dist2(b,c) > 10.0)  ) {
            printf("suspicious triangle %d (%f,%f)-(%f,%f)-(%f,%f)\n",i, a[0],a[1],b[0],b[1],c[0],c[1]);
        }
        if ( (dist2(a,c) == 0.0) ||(dist2(b,a) == 0.0) || (dist2(b,c) == 0.0)  ) {
            printf("suspicious triangle %d (%f,%f)-(%f,%f)-(%f,%f)\n",i, a[0],a[1],b[0],b[1],c[0],c[1]);
        }
        */

        //float dot00 = dot2(v0, v0);
       // float dot01 = dot2(v0, v1);
       // float dot02 = dot2(v0, v2);
       // float dot11 = dot2(v1, v1);
       // float dot12 = dot2(v1, v2);

       // float invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01) ;
       if (( v0[0]*v1[1] - v0[1]*v1[0] )!=0.0) {
        float invDenom = 1.0 / ( v0[0]*v1[1] - v0[1]*v1[0] ) ;

      //  float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
      //  float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
      //  float w = 1.0 - u - v;

        float w = (v2[0] * v1[1] - v2[1]*v1[0]) * invDenom;
        float v = (v0[0] * v2[1] - v0[1]*v2[0]) * invDenom;
        float u = 1.0 - v - w;


        if ((u>=0.0) && (v>=0.0) && (w>=0.0)) { 
            uv[0] =      u*x[0] + v*y[0] + w*z[0];
            uv[1] =1.0*( u*x[1] + v*y[1] + w*z[1]);
#ifdef DEBUG_MI
            printf("hit %f %f <- %f %f ; %d %f %f %f\n", uv[0],uv[1],savep[0],savep[1],i,u,v,w);
#endif
            *hit_count += 1;
            return;
        }

        float mm = max(fabs(0.5-u),fabs(0.5-v),fabs(0.5-w));
        if (mm<best_mm) {
            best_mm = mm;
            best_u = u;
            best_v = v;
            best_w = w;
            best_triangle = i;
        }
        }
        it = next_index(it);
    }

    i = best_triangle;
    if ((i==-1)) { 
        uv[0] = 0.0; // 1.0*(     u*x[0] + v*y[0] + w*z[0]);
        uv[1] = 0.0; // 1.0*( u*x[1] + v*y[1] + w*z[1]);
        return;
    } 

    float* a = & mt->points[ mt->triangles[i*3]*2    ];
    float* b = & mt->points[ mt->triangles[i*3+1]*2 ];
    float* c = & mt->points[ mt->triangles[i*3+2]*2 ];

    float u = best_u;
    float v = best_v;
    float w = best_w;

    float* x = & mt->tpoints[ mt->ttriangles[i*3]*2   ];
    float* y = & mt->tpoints[ mt->ttriangles[i*3+1]*2 ];
    float* z = & mt->tpoints[ mt->ttriangles[i*3+2]*2 ];

    uv[0] = 1.0*(     u*x[0] + v*y[0] + w*z[0]);
    uv[1] =1.0*( u*x[1] + v*y[1] + w*z[1]);


#ifdef DEBUG_MI

    printf("best mm: %f %f ; %d %f %f %f ( %f , %f ) %f %f ; %f %f ; %f %f \n",uv[0], uv[1], i,u,v,w, savep[0],savep[1], a[0], a[1], b[0], b[1], c[0], c[1]);
#endif

    return;
    

}


#define DEBUG_CLIPPING 0

void clip_triangle( float* a, float* b, float* c, float* d, Plane* p, int* nt, int* altered ) {
    *altered=0;
	float da = signed_distance_to_plane(a,p);
	float db = signed_distance_to_plane(b,p);
	float dc = signed_distance_to_plane(c,p);
	if (da>=0.0 && db>=0.0 && dc>=0.0) { *nt = 1; return; }
	if (da<=0.0 && db<=0.0 && dc<=0.0) { *nt = 0; return; }

#if DEBUG_CLIPPING
    printf("clip plane %f %f %f ; %f\n",p->n[0], p->n[1], p->n[2], p->d );
    printf("da = %f\n",da );
    printf("db = %f\n",db );
    printf("dc = %f\n",dc );
#endif

	if (da>0.0 && db<=0.0 && dc<=0.0 ) {   // only a survives;  
		segment_plane_intersection(a,b,b,p);  // replace b with ab-plane intersection 
		segment_plane_intersection(a,c,c,p);  // replace c with ac-plane intersection 
		*nt = 1; 
        *altered=1;
		return; 
	}

	if (db>0.0 && da<=0.0 && dc<=0.0 ) {   // only b survives;  
		segment_plane_intersection(a,b,a,p);  // replace a with ab-plane intersection 
		segment_plane_intersection(b,c,c,p);  // replace c with bc-plane intersection 
		*nt = 1; 
        *altered=1;
		return; 
	}

	if (dc>0.0 && da<=0.0 && db<=0.0 ) {   // only c survives;  
		segment_plane_intersection(a,c,a,p);  // replace a with ab-plane intersection 
		segment_plane_intersection(b,c,b,p);  // replace b with bc-plane intersection 
		*nt = 1; 
        *altered=1;

		return; 
	}


	if (da<0.0 && db>=0.0 && dc>=0.0 ) {   // only a rejected; 
    //printf("a rejected\n"); 
        float a0[3];
        set(&a0[0],a);
		segment_plane_intersection(a,b,a,p);  // replace b with ab-plane intersection 

		segment_plane_intersection(a0,c,d,p);  // set the new vertex d 
		set(&d[3],a);
		set(&d[6],c);

		*nt = 2; 
        *altered=1;

		return; 
	}


	if (dc<0.0 && da>=0.0 && db>=0.0 ) {   // only c rejected;  
    //    printf("c rejected\n"); 

        float c0[3];
        set(c0,c);
//        printf("copy of c= %.2f %.2f %.2f\n",c0[0],c0[1],c0[2]);
		segment_plane_intersection(c,b,c,p);  // replace c with cb-plane intersection 

//        printf("new vertex1 %.2f %.2f %.2f between (%.2f %.2f %.2f) and (%.2f %.2f %.2f)\n", c[0], c[1], c[2],
//        c0[0], c0[1], c0[2],b[0], b[1], b[2]
 //       );

		segment_plane_intersection(c0,a,d,p);  // set the new vertex d  
   //     printf("new vertex2 %.2f %.2f %.2f between (%.2f %.2f %.2f) and (%.2f %.2f %.2f)\n", d[0], d[1], d[2],
   //     c0[0], c0[1], c0[2],a[0], a[1], a[2]
    //    );
		set(&d[3],a);
		set(&d[6],c);

		*nt = 2; 
        *altered=1;

		return; 
	}

	if (db<0.0 && da>=0.0 && dc>=0.0 ) {   // only b rejected;  
    //    printf("b rejected\n"); 

        float b0[3];
        set(&b0[0],b);

		segment_plane_intersection(a,b,b,p);  // replace c with cb-plane intersection 

		segment_plane_intersection(b0,c,d,p);  // set the new vertex d  
		set(&d[3],c);
		set(&d[6],b);

		*nt = 2; 
        *altered=1;

		return; 	
	}

}


bool triangle_spans_plane( float* a, float* b, float* c, Plane* p ) {
	if ( segment_spans_plane(a,b,p) ) {return true;}
	if ( segment_spans_plane(b,c,p) ) {return true;}
	if ( segment_spans_plane(c,a,p) ) {return true;}
	if ((signed_distance_to_plane(a,p) > 0.0) && (signed_distance_to_plane(b,p) < 0.0)) {return true;}
	return false;
}

void segment_plane_intersection( float* a, float* b, float* q, Plane* p ) {
    

//    printf("find segment intersection on (%.2f %.2f %.2f) - (%.2f %.2f %.2f)\n", 
//        a[0], a[1], a[2],b[0], b[1], b[2]
//    );


	float v[3];
	subtract(v,b,a);  // v = b-a

//    printf("difference vector = (%.2f %.2f %.2f) \n", 
//       v[0], v[1], v[2]
//    );

	float t = (-p->d - dot( a , &p->n[0] )) / dot( v, &p->n[0] ) ;
  //  printf("t=%f\n",t);
	q[0] = a[0] + t*v[0];
	q[1] = a[1] + t*v[1];
	q[2] = a[2] + t*v[2];

    //    printf("result = (%.2f %.2f %.2f) \n",  q[0],q[1],q[2]);


//    printf("find segment intersection on (%.2f %.2f %.2f) - (%.2f %.2f %.2f)\n", 
//        a[0], a[1], a[2],b[0], b[1], b[2]
//    );

}


bool segment_spans_plane( float* a, float* b, Plane* p ) {
	if ((signed_distance_to_plane(a,p) < 0.0) && (signed_distance_to_plane(b,p) > 0.0)) {return true;}
	if ((signed_distance_to_plane(a,p) > 0.0) && (signed_distance_to_plane(b,p) < 0.0)) {return true;}
	return false;
}

float signed_distance_to_plane(float* a, Plane* p) {
	return (dot(a,&(p->n[0]))+ p->d);
}


int add_side_quad( float* triangle ,  Plane* p, float* new_triangles, float droplevel ) {
   //printf("add sidequad\n"); fflush(stdout);
    float *a = &triangle[0];
    float *b = &triangle[3];
    float *c = &triangle[6];

    float *x;
    float *y;

    float z[3];
    float w[3];

    int n_on_plane=0;
    float da = signed_distance_to_plane(a,p);  if (da==0.0) {n_on_plane++;}
    float db = signed_distance_to_plane(b,p);  if (db==0.0) {n_on_plane++;}
    float dc = signed_distance_to_plane(c,p);  if (dc==0.0) {n_on_plane++;}

    if (n_on_plane == 0) {return 0;}
    if (n_on_plane == 1) {return 0;}
    if (n_on_plane == 3) {printf("warning: new triangle in the plane?"); return 0;}

    if (da==0.0 && db==0.0) {  x = a;  y = b;    }
    if (db==0.0 && dc==0.0) {  x = b;  y = c;    }
    if (dc==0.0 && da==0.0) {  x = c;  y = a;    }

    set(&z[0],x); z[2]= droplevel;
    set(&w[0],y); w[2] = droplevel;
 
    if (1) {
    set( &new_triangles[0] , x );
    set( &new_triangles[3] , z );
    set( &new_triangles[6] , y );

    set( &new_triangles[9] , y );
    set( &new_triangles[12] , z );
    set( &new_triangles[15] , w );
    } else {
    set( &new_triangles[18] , x );
    set( &new_triangles[21] , y );
    set( &new_triangles[24] , z );

    set( &new_triangles[27] , z );
    set( &new_triangles[30] , y );
    set( &new_triangles[33] , w );
    }
    return(2);
}
