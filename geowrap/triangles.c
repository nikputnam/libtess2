
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "triangles.h"


MeshTriangles* parse_triangles(char* filename, float width) {

    MeshTriangles* mt = malloc(sizeof(MeshTriangles));

    FILE* file = fopen(filename, "r"); 
    char line[1024];

    printf("%s\n",filename);
    //printf("%x\n",file);

    int npoints = 0;
    int ntriangles = 0;
    int ntpoints = 0;
    while (fgets(line, sizeof(line), file)) {
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
        if (strncmp(line,"vn ",3)==0) { 
            //printf("%s", line); 
            ntpoints++;
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

    // x3 because we're replicating the interpolation mesh 
    mt->npoints = 3*npoints;
    mt->ntpoints = 3*ntpoints;
    mt->ntriangles = 3*ntriangles;
    mt->npaths = 3*ntriangles;

    printf("npoints: %d\n",npoints);
    printf("ntpoints: %d\n",ntpoints);
    printf("ntriangles: %d\n",ntriangles);

    fclose(file);

    mt->points = malloc( sizeof(float)*npoints*2 *3 );  assert(mt->points != 0);
    mt->tpoints = malloc( sizeof(float)*ntpoints*2 *3 );assert(mt->tpoints != 0);
    mt->triangles = malloc( sizeof(int)*ntriangles*3 *3 );assert(mt->triangles != 0);
    mt->ttriangles = malloc( sizeof(int)*ntriangles*3 *3 );assert(mt->ttriangles != 0);
    mt->paths = malloc( sizeof(float)*ntriangles*6 *3 );assert(mt->paths != 0);

    file = fopen(filename, "r");   assert(file);
    printf("x\n"); fflush(stdout);
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    int m=0;
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
            int a,b,c,x,y,z;
            
            sscanf(line,"f %d/%d %d/%d %d/%d",&a,&x,&b,&y,&c,&z);
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

        k=0;
        for (int i=0; i<mt->ntriangles;i++ ) {
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3]*2  ];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3]*2+1];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+1]*2 ];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+1]*2 +1];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+2]*2 ];
                mt->paths[k++] = mt->tpoints[ mt->triangles[i*3+2]*2 +1];
        }

        assert(k==ntriangles*6*3);

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

#define DEBUG_MI 1

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
