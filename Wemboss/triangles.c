
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "triangles.h"
#include "vectorops.h"



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


void apply_interpolation_to_mesh( MeshTriangles* t, MeshTriangles* interpolater, float zscale ) {

    for (int i=0;i<t->nxpoints;i++) {

       // printf("%f %f %f\n", t->xpoints[3*i], t->xpoints[3*i+1], t->xpoints[3*i+2]);
        mesh_interpolation(interpolater, & t->xpoints[3*i], & t->xpoints[3*i] );
        t->xpoints[3*i+2] = zscale*t->xpoints[3*i+2]; 
        //t->xpoints[3*i+1] = t->xpoints[3*i+1] * -1.0;

        //printf("                                   -> %f %f %f\n", t->xpoints[3*i], t->xpoints[3*i+1], t->xpoints[3*i+2]);

        if (i%10000 == 0) { printf("%d / %d = %f\n",i,t->nxpoints, ((float)i) / t->nxpoints) ; }
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

MeshTriangles* parse_triangles_internal(char* filename, float width, int with_normals, int reduplicate) {

    MeshTriangles* mt = malloc(sizeof(MeshTriangles));

    FILE* file = fopen(filename, "r"); 
    char line[1024];

    printf("%s\n",filename);
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
            //printf("%s", line); 
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

    fclose(file);

    mt->points = malloc( sizeof(float)*(mt->npoints)*2  );  assert(mt->points != 0);
    mt->tpoints = malloc( sizeof(float)*(mt->ntpoints)*2  );assert(mt->tpoints != 0);
    mt->xpoints = malloc( sizeof(float)*(mt->nxpoints)*3 );assert(mt->xpoints != 0);
    mt->triangles = malloc( sizeof(int)*(mt->ntriangles)*3  );assert(mt->triangles != 0);
    mt->ttriangles = malloc( sizeof(int)*(mt->ntriangles)*3  );assert(mt->ttriangles != 0);
  //  mt->xtriangles = malloc( sizeof(int)*ntriangles*3 );assert(mt->xtriangles != 0);
    mt->paths = malloc( sizeof(float)*ntriangles*6 );assert(mt->paths != 0);

    file = fopen(filename, "r");   assert(file);
    printf("x\n"); fflush(stdout);
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    int m=0;
    int xi=0;
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
            
            if (with_normals) {
                sscanf(line,"f %d/%d/%d %d/%d/%d %d/%d/%d",&a,&x,&n1,&b,&y,&n2,&c,&z,&n3);
            } else {
                sscanf(line,"f %d/%d %d/%d %d/%d",&a,&x,&b,&y,&c,&z);
            }
            //flip triangles
//            sscanf(line,"f %d/%d %d/%d %d/%d",&c,&z,&b,&y,&a,&x);
//            sscanf(line,"f %d/%d %d/%d %d/%d",&c,&x,&b,&y,&a,&z);
            //printf("%d %d %d %s\n",a,b,c,line);
            mt->triangles[j++]=a-1;
            mt->triangles[j++]=b-1;
            mt->triangles[j++]=c-1;

            mt->ttriangles[m++]=x-1;
            mt->ttriangles[m++]=y-1;
            mt->ttriangles[m++]=z-1;

/*
            mt->paths[k++] = mt->tpoints[2*(x-1)];
            mt->paths[k++] = mt->tpoints[2*(x-1)+1];
            mt->paths[k++] = mt->tpoints[2*(y-1)];
            mt->paths[k++] = mt->tpoints[2*(y-1)+1];
            mt->paths[k++] = mt->tpoints[2*(z-1)];
            mt->paths[k++] = mt->tpoints[2*(z-1)+1];
*/
//            mt->points[i++]=y;
 //           ntriangles++;
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
        k=0;
        for (int i=0; i<mt->ntriangles;i++ ) {
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3]*2  ];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3]*2+1];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+1]*2 ];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+1]*2 +1];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+2]*2 ];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+2]*2 +1];
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


float dot2(float* x, float* y){
    return (x[0]*y[0] + x[1]*y[1]);
}

void diff(float* x, float* y, float* z) {
    z[0] = x[0]-y[0];
    z[1] = x[1]-y[1];
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

void mesh_interpolation(MeshTriangles* mt, float* p, float* uv) {

//    uv[0]=0.0;
//    uv[1]=0.0;
//    return;

    float best_mm = 1000000000.0;
    float best_u = 10.0;
    float best_v = 10.0;
    float best_w = 10.0;
    int best_triangle = 0;

#ifdef DEBUG_MI
    float savep[2] = {p[0], p[1]};
#endif

    for (int i=0; i<mt->ntriangles;i++ ) {
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

        //float dot00 = dot2(v0, v0);
       // float dot01 = dot2(v0, v1);
       // float dot02 = dot2(v0, v2);
       // float dot11 = dot2(v1, v1);
       // float dot12 = dot2(v1, v2);

       // float invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01) ;
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

    int i = best_triangle;

        float* a = & mt->points[ mt->triangles[i*3]*2    ];
        float* b = & mt->points[ mt->triangles[i*3+1]*2 ];
        float* c = & mt->points[ mt->triangles[i*3+2]*2 ];

    float u = best_u;
    float v = best_v;
    float w = best_w;

        float* x = & mt->tpoints[ mt->ttriangles[i*3]*2    ];
        float* y = & mt->tpoints[ mt->ttriangles[i*3+1]*2 ];
        float* z = & mt->tpoints[ mt->ttriangles[i*3+2]*2 ];

        uv[0] = 1.0*(     u*x[0] + v*y[0] + w*z[0]);
        uv[1] =1.0*( u*x[1] + v*y[1] + w*z[1]);

#ifdef DEBUG_MI

    printf("best mm: %f %f ; %d %f %f %f ( %f , %f ) %f %f ; %f %f ; %f %f \n",uv[0], uv[1], i,u,v,w, savep[0],savep[1], a[0], a[1], b[0], b[1], c[0], c[1]);
#endif

    return;
    

}
