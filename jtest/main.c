#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "frozen.h"
#include "surface.h"
#include "catmullrom.h"


int print_array(struct json_out * o, va_list *ap) {

  int n = va_arg(*ap, int);
  float* xy = va_arg(*ap, float*);
  printf("print array got %d\n",n);
  for (int i=0; i<n-1 ; i++) {
    printf("xy: %f %f\n",xy[0], xy[1]);
    json_printf(o, "{x: %f, y:%f},", xy[0],xy[1]);
    xy++; xy++;
  }
  json_printf(o, "{x: %f, y:%f}", xy[0],xy[1]);

  return 0;
}


int main(int argc, char** argv) {
  //  printf("hello world\n");
    pottoy_spec_t spec;  

    read_spec(argv[1], &spec );
    float scaled_height = atof(argv[3]);

    //read_spec(argv[1], &spec );
    //printf("read spec with %d points\n",spec.npoints);

 		float bbox[6];
		surface_obj_bbox(&bbox[0], &spec) ;
		printf("bbox: x:  %f - %f\n", bbox[0], bbox[3]);
		printf("bbox: y:  %f - %f\n", bbox[1], bbox[4]);
		printf("bbox: z:  %f - %f\n", bbox[2], bbox[5]);

    for (int i = 0; i< spec.npoints*2; i+=2) {
      spec.points[i+1] -= bbox[1];

      spec.points[i] *= scaled_height/(bbox[4]-bbox[1]);
      spec.points[i+1] *= scaled_height/(bbox[4]-bbox[1]);
      printf("xy %f %f\n",spec.points[i], spec.points[i+1]);
    }


//    for (int i = 0; i< spec.npoints*2; i+=2) {
 //       printf("xy %f %f\n",spec.points[i], spec.points[i+1]);
 //   }

    printf("squash %f\n",spec.squash);
    printf("twist %f\n",spec.twist);

    float x[2];
    for (float t=0; t<1.0; t+=0.01) {
        catmullrom2(t,&spec.points[0], spec.npoints, &x[0]);
        //printf("%f %f %f\n", t,x[0],x[1]);
    }

  int n_sectors = 90;
  if (spec.facet) {
    int spf = n_sectors / spec.n_facets ;
    if (spf < 10 ) {spf = 10;}
    n_sectors = spec.n_facets * (spf + 1);
    printf("n_sectors = %d\n", n_sectors );
  }

  write_surface_obj0(argv[2], &spec, n_sectors, 50);     



  json_fprintf("out.json","{ n_facets:%d, facet:%B, npoints:%d, puff:%f, squash:%f, twist:%f, points:[%M] }", 
      spec.n_facets, spec.facet, spec.npoints, spec.puff, spec.squash, spec.twist, print_array, spec.npoints, spec.points
  );


}

