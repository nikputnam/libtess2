typedef struct MeshTriangles
{
	int ntriangles;				
	int npoints;				
	int ntpoints;				
	float* points;			
	float* tpoints;			
	int* triangles;		
	int* ttriangles;		
    float* paths;	
} MeshTriangles;


MeshTriangles* parse_triangles(char* filename, float width) ;
void mesh_interpolation(MeshTriangles* mt, float* xy, float* uv) ;
