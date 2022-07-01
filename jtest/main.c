#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "frozen.h"
#include "surface.h"
#include "catmullrom.h"

int main(int argc, char** argv) {
  //  printf("hello world\n");
    pottoy_spec_t spec;  

    read_spec(argv[1], &spec );
    //printf("read spec with %d points\n",spec.npoints);

 
    for (int i = 0; i< spec.npoints*2; i+=2) {
        printf("xy %f %f\n",spec.points[i], spec.points[i+1]);
    }

    printf("squash %f\n",spec.squash);
    printf("twist %f\n",spec.twist);

    float x[2];
    for (float t=0; t<1.0; t+=0.01) {
        catmullrom2(t,&spec.points[0], spec.npoints, &x[0]);
        printf("%f %f %f\n", t,x[0],x[1]);
    }

    write_surface_obj2(argv[2], &spec);     

  //  free(buffer);

}

