
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
        if (strncmp(line,"v ",2)==0) { 
            //printf("%s", line); 
            npoints++;
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

    mt->npoints = npoints;
    mt->ntpoints = ntpoints;
    mt->ntriangles = ntriangles;

    printf("npoints: %d\n",npoints);
    printf("ntriangles: %d\n",ntriangles);
    fclose(file);

    mt->points = malloc( sizeof(float)*npoints*2 );
    mt->tpoints = malloc( sizeof(float)*ntpoints*2 );
    mt->triangles = malloc( sizeof(int)*ntriangles*3 );
    mt->ttriangles = malloc( sizeof(int)*ntriangles*3 );
    mt->paths = malloc( sizeof(float)*ntriangles*6 );

    file = fopen(filename, "r"); 

    int i=0;
    int j=0;
    int k=0;
    int l=0;
    int m=0;
    while (fgets(line, sizeof(line), file)) {
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
        if (strncmp(line,"v ",2)==0) { 
            float x;
            float y;
            //printf("%s", line); 
            sscanf(line,"v %f %f",&x,&y);
            //printf("%f %f %s\n",x,y,line);
            mt->points[i++]=x*width;
            mt->points[i++]=y*width;
//            mt->points[i++]=x;
//            mt->points[i++]=y;
        }

        if (strncmp(line,"vt ",3)==0) { 
            float x;
            float y;
            //printf("%s", line); 
            sscanf(line,"vt %f %f",&x,&y);
            //printf("%f %f %s\n",x,y,line);
            mt->tpoints[l++]=x;
            mt->tpoints[l++]=y;
        }

        if (strncmp(line,"f ",2)==0) { 
            int a,b,c,x,y,z;
            
            //sscanf(line,"f %d/%d %d/%d %d/%d",&a,&x,&b,&y,&c,&z);
            //flip triangles
            sscanf(line,"f %d/%d %d/%d %d/%d",&c,&z,&b,&y,&a,&x);
            //printf("%d %d %d %s\n",a,b,c,line);
            mt->triangles[j++]=a-1;
            mt->triangles[j++]=b-1;
            mt->triangles[j++]=c-1;

            mt->ttriangles[m++]=x-1;
            mt->ttriangles[m++]=y-1;
            mt->ttriangles[m++]=z-1;

            mt->paths[k++] = mt->tpoints[2*(x-1)];
            mt->paths[k++] = mt->tpoints[2*(x-1)+1];
            mt->paths[k++] = mt->tpoints[2*(y-1)];
            mt->paths[k++] = mt->tpoints[2*(y-1)+1];
            mt->paths[k++] = mt->tpoints[2*(z-1)];
            mt->paths[k++] = mt->tpoints[2*(z-1)+1];

//            mt->points[i++]=y;
//            printf("%s", line); 
 //           ntriangles++;
        }
    }
    printf("%d %d\n",i,j);

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


void mesh_interpolation(MeshTriangles* mt, float* p, float* uv) {

//    uv[0]=0.0;
//    uv[1]=0.0;
//    return;

    float best_mm = 1000000000.0;
    float best_u = 10.0;
    float best_v = 10.0;
    float best_w = 10.0;
    int best_triangle = 0;

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
            uv[0] = u*x[0] + v*y[0] + w*z[0];
            uv[1] = u*x[1] + v*y[1] + w*z[1];
            printf("hit! %f %f <- %f %f ; %d %f %f %f\n", uv[0],uv[1],p[0],p[1],i,u,v,w);
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

        uv[0] = u*x[0] + v*y[0] + w*z[0];
        uv[1] = u*x[1] + v*y[1] + w*z[1];

    printf("best mm: %f %f ; %d %f %f %f ( %f , %f ) %f %f ; %f %f ; %f %f \n",uv[0], uv[1], i,u,v,w, p[0],p[1], a[0], a[1], b[0], b[1], c[0], c[1]);


    return;
    

}
