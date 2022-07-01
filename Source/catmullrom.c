#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vectorops.h"
#include "catmullrom.h"


void catmullrom2(float t, float* points, int npoints, float* x) {

    float p0[2];
    float p1[2];
    float p2[2];
    float p3[2];
    
    int i = (int) (t*((float) (npoints-1)));
    if (i==(npoints-1)) {i--;}  // in case t==1.0000;

    float tp = (t-(1.0/(npoints-1))*((float) i))*(npoints-1) ;
    float alpha = 0.5;    // 0.5 for "centripital"

  //  printf("catmull rom debug t,i,tp npoints: %f %d %f %d\n", t, i, tp, npoints);

    if (i==0) {
        subtract2(&p0[0], &points[0], &points[2]);  //  v0 - v1
        add2(&p0[0], &p0[0], &points[0]);
    } else {
        set2(&p0[0], &points[(i-1)*2] );
    }

    set2(&p1[0], &points[i*2] );
    set2(&p2[0], &points[(i+1)*2] );

    if (i==(npoints-2)) { // the final segment
        subtract2(&p3[0], &points[(i+1)*2], &points[i*2]);  //  v0 - v1
        add2(&p3[0], &p3[0], &points[(i+1)*2]);
    } else {
        set2(&p3[0], &points[(i+2)*2] );
    }

  //  printf("debug p0 %f %f %f\n", p0[0], p0[1], t);
  //  printf("debug p1 %f %f %f\n", p1[0], p1[1], t);
  //  printf("debug p2 %f %f %f\n", p2[0], p2[1], t);
  //  printf("debug p3 %f %f %f\n", p3[0], p3[1], t);

    catmullrom2_segment(tp,p0,p1,p2,p3,alpha,x);

}

void catmullrom2_segment(float tp, float* p0, float* p1, float* p2, float* p3, float alpha, float* x) {

    float A1[2];
    float A2[2];
    float A3[2];

    float B1[2];
    float B2[2];

    float t0 = 0.0;
    float t1 = t0 + pow(dist2(p0,p1),alpha);
    float t2 = t1 + pow(dist2(p1,p2),alpha);
    float t3 = t2 + pow(dist2(p2,p3),alpha);

    float t = t1 + tp*(t2-t1);

    weighted_sum2( &A1[0], ((t1-t)/(t1-t0)), ((t-t0)/(t1-t0)), p0, p1  );
    weighted_sum2( &A2[0], ((t2-t)/(t2-t1)), ((t-t1)/(t2-t1)), p1, p2  );
    weighted_sum2( &A3[0], ((t3-t)/(t3-t2)), ((t-t2)/(t3-t2)), p2, p3  );

    weighted_sum2( &B1[0], ((t2-t)/(t2-t0)), ((t-t0)/(t2-t0)), &A1[0], &A2[0]  );
    weighted_sum2( &B2[0], ((t3-t)/(t3-t1)), ((t-t1)/(t3-t1)), &A2[0], &A3[0]  );

    weighted_sum2( x     , ((t2-t)/(t2-t1)), ((t-t1)/(t2-t1)), &B1[0], &B2[0]  );

}