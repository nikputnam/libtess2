#include <stdlib.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "surface.h"

#include "frozen.h"
#include "catmullrom.h"


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
    catmullrom2(uv[0],&spec->points[0], spec->npoints, &point[0]);

	float r = point[0];
	float y = point[1];

    xyz[0] = r*cos( uv[1]* M_PI *2.0) ;
    xyz[1] = y;
    xyz[2] = spec->squash*r*sin( uv[1] * M_PI *2.0);

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
    float v0 = uv[1];
    float pufff =   (1.0-exp(-u*20.0))*(1-exp(-1.0*(1.0-u)*20.0)); //min(1,10*u,10*(1-u)) ;
	float v = v0 + u * (spec->twist) ;

    double ip;
	float phase = modf( ( (v0) * ((float) spec->n_facets)),&ip);
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


// output an obj file augmented with some vertices having fixed UV coordiantes based on the parametric fn X()
int write_surface_obj2(char* fn, pottoy_spec_t* spec) {

    FILE* objf = fopen(fn, "wt");
    float xyz[3];
    float uv[2];

    int n_sectors=90;
    int n_levels=50;

    for(int j=0;j<n_sectors;j++) {
        double v = (((double) j)/((double) (n_sectors-1)));
        printf("v: %f\n",v);
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
    
    /* reset the file position indicator to 
    the beginning of the file */
    fseek(infile, 0L, SEEK_SET);	
    
    /* grab sufficient memory for the 
    buffer to hold the text */
    buffer = (char*)calloc(numbytes, sizeof(char));	
    
    /* memory error */
    if(buffer == NULL)
        return 1;
    
    /* copy all the text into the buffer */
    fread(buffer, sizeof(char), numbytes, infile);
    fclose(infile);
    


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

//    printf("n_facets: %d\n",a);
 //   printf("facet: %d\n",d);

    spec->points = (float*)calloc( spec->npoints, 2*sizeof(float) );


      json_scanf(buffer, strlen(buffer), "{ points:%M  }",
             scan_array, spec->points );

/*    for (int i = 0; i< n_points*2; i++) {
        printf("xy %f\n",points[i]);
    }
*/
    free(buffer);
    return(0);

}
