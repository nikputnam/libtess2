#include <stdlib.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "surface.h"

void X(float* xyz, float* vtheta) {
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
}
