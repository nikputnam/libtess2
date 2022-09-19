typedef struct MeshTriangles
{
	int ntriangles;				
	int npaths;				
	int npoints;				
	int nxpoints;				
	int ntpoints;				
	float* points;			
	float* xpoints;			
	float* tpoints;			
	int* triangles;		
	int* ttriangles;		
    float* paths;	
} MeshTriangles;


MeshTriangles* parse_triangles(char* filename, float width) ;
MeshTriangles* parse_triangles_internal(char* filename, float width,int with_normals) ;
MeshTriangles* parse_triangles_with_normals(char* filename, float width) ;
void mesh_interpolation(MeshTriangles* mt, float* xy, float* uv) ;
void write_to_stl( MeshTriangles* t, FILE* stlfile ) ;


void print_triangle_raw( float* v1, float* v2, float* v3 , FILE* fp) ;
void triangle_normal(float* v1, float* v2, float* v3, float *n) ;

void apply_interpolation_to_mesh( MeshTriangles* texture, MeshTriangles* mt , float zscale);
void apply_transform_to_mesh(MeshTriangles*  texture,  void(*trnsfrm)(float*, float*) );
