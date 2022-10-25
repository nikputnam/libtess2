#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
//#include "nanosvg.h"
//#include "bez.h"
//#include "tess.h"
//#include "tesselator.h"
#include "vectorops.h"
#include "triangles.h"

int main(int argc, char *argv[]) {

		if (argc!=3) {
			printf("usage: cylzip <in.obj> <out.obj>\n");
			exit(0);
		}

		Plane p = (Plane){ {0,1,0}, 0.00 };

        MeshTriangles* mesh = parse_triangles_with_normals(argv[1],(float) 1.0);

		int* flags= calloc(mesh->nxpoints,sizeof(int));
//		memset(flags,0,mesh->nxpoints);

		printf("parsed triangles\n");

		for (int i=0; i<mesh->ntriangles; i++) {
			int v1 = mesh->triangles[i*3];
			int v2 = mesh->triangles[i*3+1];
			int v3 = mesh->triangles[i*3+2];

			float* x1 = &mesh->xpoints[v1*3];
			float* x2 = &mesh->xpoints[v2*3];
			float* x3 = &mesh->xpoints[v3*3];

			float d1 =  signed_distance_to_plane(x1, &p);
			float d2 =  signed_distance_to_plane(x2, &p);
			float d3 =  signed_distance_to_plane(x3, &p);

			int n_lz = 0;
			int n_gez = 0;
			if (d1<0.0) {n_lz++;} 
			if (d2<0.0) {n_lz++;} 
			if (d3<0.0) {n_lz++;} 
			if (d1>=0.0) {n_gez++;} 
			if (d2>=0.0) {n_gez++;} 
			if (d3>=0.0) {n_gez++;} 
			printf("zz %f %f %f\n",d1,d2,d3);
			if ((n_lz>0 && n_gez>0 && x1[0]>0.0))  { printf("#################"); }
			if ((n_lz>0 && n_gez>0 && x1[0]>0.0))  {
				printf("f %d %d %d\n",  mesh->triangles[i*3]+1,  mesh->triangles[i*3+1]+1,  mesh->triangles[i*3+2]+1  );
				if (d1>=0) { flags[v1]=1; printf("v %f %f %f  %d\n", mesh->xpoints[v1*3],  mesh->xpoints[v1*3+1],  mesh->xpoints[v1*3+2], v1 ); }
				if (d2>=0) { flags[v2]=1; printf("v %f %f %f  %d\n", mesh->xpoints[v2*3],  mesh->xpoints[v2*3+1],  mesh->xpoints[v2*3+2], v2 ); }
				if (d3>=0) { flags[v3]=1; printf("v %f %f %f  %d\n", mesh->xpoints[v3*3],  mesh->xpoints[v3*3+1],  mesh->xpoints[v3*3+2], v3 ); }
				printf( "d1 %f\n",  signed_distance_to_plane(x1, &p) );
				printf( "d2 %f\n",  signed_distance_to_plane(x2, &p) );
				printf( "d3 %f\n",  signed_distance_to_plane(x3, &p) );
			}
		}

		FILE* out = fopen(argv[2],"wt");

		int n = mesh->nxpoints;
		for (int i=0;i<mesh->nxpoints;i++) {
			if (flags[i]) {
				flags[i] = n++;
			}
		}

		for (int i=0;i<mesh->nxpoints;i++) {
			if (flags[i]) {
				if (flags[i]==mesh->nxpoints) {
					fprintf(out,"v %f %f %f fix 0.0 0.0\n",mesh->xpoints[i*3], mesh->xpoints[i*3+1], mesh->xpoints[i*3+2] );
				} else {
					fprintf(out,"v %f %f %f zip %d\n",mesh->xpoints[i*3], mesh->xpoints[i*3+1], mesh->xpoints[i*3+2], flags[i]-mesh->nxpoints );
				}
			} else {
				fprintf(out,"v %f %f %f\n",mesh->xpoints[i*3], mesh->xpoints[i*3+1], mesh->xpoints[i*3+2]);
			}
		}

		for (int i=0;i<mesh->nxpoints;i++) {
			if (flags[i]) {
				if (flags[i]==mesh->nxpoints) {
					fprintf(out,"v %f %f %f fix 1.0 0.0\n",mesh->xpoints[i*3], mesh->xpoints[i*3+1], mesh->xpoints[i*3+2] );
				} else {
					fprintf(out,"v %f %f %f zip %d\n",mesh->xpoints[i*3], mesh->xpoints[i*3+1], mesh->xpoints[i*3+2], flags[i]-mesh->nxpoints );
				}
			}
		}


		for (int i=0;i<mesh->nxpoints;i++) {
			if (flags[i]) {
				printf("zip: %d -> %d\n",i,flags[i]);
			}
		}


		for (int i=0; i<mesh->ntriangles; i++) {
			int v1 = mesh->triangles[i*3];
			int v2 = mesh->triangles[i*3+1];
			int v3 = mesh->triangles[i*3+2];

			float* x1 = &mesh->xpoints[v1*3];
			float* x2 = &mesh->xpoints[v2*3];
			float* x3 = &mesh->xpoints[v3*3];

			float d1 =  signed_distance_to_plane(x1, &p);
			float d2 =  signed_distance_to_plane(x2, &p);
			float d3 =  signed_distance_to_plane(x3, &p);

			int n_lz = 0;
			int n_gez = 0;
			if (d1<0.0) {n_lz++;} 
			if (d2<0.0) {n_lz++;} 
			if (d3<0.0) {n_lz++;} 
			if (d1>=0.0) {n_gez++;} 
			if (d2>=0.0) {n_gez++;} 
			if (d3>=0.0) {n_gez++;} 
			//printf("zz %f %f %f\n",d1,d2,d3);
			//if ((n_lz>0 && n_gez>0 && x1[0]>0.0))  { printf("#################"); }
			if ((n_lz>0 && n_gez>0 && x1[0]>0.0))  {
				printf("f %d %d %d xxxx\n",  mesh->triangles[i*3]+1,  mesh->triangles[i*3+1]+1,  mesh->triangles[i*3+2]+1  );
				if (d1>=0) { v1 = flags[v1]; }
				if (d2>=0) { v2 = flags[v2]; }
				if (d3>=0) { v3 = flags[v3]; }
			//	fprintf(out, "f %d %d %d zip\n", v1+1, v2+1, v3+1 );
			}
			fprintf(out, "f %d %d %d\n", v1+1, v2+1, v3+1 );
		}


		fclose(out);

		printf("done\n");

}