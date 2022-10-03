

#ifndef _meshindex2d_h_
#define _meshindex2d_h_

typedef struct meshindex
{
    float* points;
    int*   triangles;
    int*   bincounts;
    int*   bins;

    float bin_width_x;
    int bin_offset_x;

    int n_xbins;
    int n_ybins;

    float bin_width_y;
    int bin_offset_y;

    int max_bin_triangles;



} meshindex;


typedef struct meshindex_it
{
    meshindex* mi;
    int a,b,i,j,status,k,t;
    int done;
} meshindex_it;


static meshindex_it* next_index(meshindex_it* it) {
    int aa = it->a + it->i ;
    int bb = it->b + it->j ;
    if (aa<0) {aa=0;}
    if (aa>=it->mi->n_xbins) {aa=it->mi->n_xbins-1;}
    if (bb<0) {bb=0;}
    if (bb>=it->mi->n_ybins) {bb=it->mi->n_xbins-1;}

    int k = it->k;
    int kmax = it->mi->bincounts[ bb * it->mi->n_xbins + aa ];
    k++;
    while (k >= kmax) {
        k=0;
        it->i += 1; if ( it->i > 1 ) { it->i = -1; it->j += 1; } 
        if (it->j > 1) { it->done = 1; return it; }


        aa = it->a + it->i ;
        bb = it->b + it->j ;

        if (aa<0) {aa=0;}
        if (aa>=it->mi->n_xbins) {aa=it->mi->n_xbins-1;}
        if (bb<0) {bb=0;}
        if (bb>=it->mi->n_ybins) {bb=it->mi->n_xbins-1;}

        kmax = it->mi->bincounts[ bb * it->mi->n_xbins + aa ];

    }
   // printf("ij %d %d (%d %d) k/kmax %d %d \n",it->i,it->j,aa,bb,k,kmax);
    it->k = k;
    if (it->j > 1) { it->done = 1; return it; }
    it->t = it->mi->bins[ (bb * it->mi->n_xbins + aa)*it->mi->max_bin_triangles + it->k ];
    return it;
}

static meshindex_it* index_iterator( float x, float y, meshindex* mi) { //}, int* k, int* r, int* c, int* status ) {

    meshindex_it* it = malloc(sizeof(meshindex_it));

    //float xx = p[0];
    //float yy = p[1];

    float xx = x;
    float yy = y;

    it->i = -1;
    it->j = -1;

//    if (status==0) {
    it->a = mi->bin_offset_x + (int) (xx/mi->bin_width_x) ;
    it->b = mi->bin_offset_y + (int) (yy/mi->bin_width_y) ;

    int aa = it->a + it->i ;
    int bb = it->b + it->j ;

    if (aa<0) {aa=0;}
    if (aa>=mi->n_xbins) {aa=mi->n_xbins-1;}
    if (bb<0) {bb=0;}
    if (bb>=mi->n_ybins) {bb=mi->n_xbins-1;}

    int k = 0 ;// mi->bincounts[ b * mi->n_xbins + a ]++;
    //mi->bins[ (b * mi->n_xbins + a)*mi->max_bin_triangles + k ] = j;

    it->mi = mi;
    it->status = 0;
    it->k = k;
    it->done = 0;
    it->t = mi->bins[ (bb * mi->n_xbins + aa)*mi->max_bin_triangles + k ];

    return it;
}

static meshindex* build_mesh_index( float* points, int* triangles, int n_triangles, int n_xbins, int n_ybins ) {

    meshindex* mi = malloc( sizeof(meshindex) );
    mi->bincounts = malloc( sizeof(int)*(1+n_xbins)*(1+n_ybins) );
    memset(mi->bincounts,0,sizeof(int)*(1+n_xbins)*(1+n_ybins));

    mi->n_xbins = n_xbins;
    mi->n_ybins = n_ybins;
    mi->points = points;
    mi->triangles = triangles;

    int max_binfill=0;

    float min_x;
    float max_x;
    float min_y;
    float max_y;

    for (int i=0;i<3;i++) {
        int a = triangles[i];

        float x = points[a*2  ];
        float y = points[a*2+1];
  //      float z = points[a*3+2];

        min_x = x;  max_x = x; 
        min_y = y;  max_y = y;   
    }

    for (int j=1;j<n_triangles;j++) {
        for (int i=0;i<3;i++) {

            int a = triangles[j*3+i];

           // printf("triangle %d, point %d, index %d\n",j,i,a);  fflush(stdout);

            float x = points[a*2  ];
            float y = points[a*2+1];
//            float z = points[a*3+2];

            min_x = min_x > x ? x : min_x; 
            max_x = max_x < x ? x : max_x; 
            min_y = min_y > y ? y : min_y; 
            max_y = max_y < y ? y : max_y; 
        }
    }

    printf("min max x %f %f\n", min_x, max_x); 
    printf("min max y %f %f\n", min_y, max_y); 


    mi->bin_width_x = ( max_x - min_x ) / n_xbins;
    mi->bin_offset_x = (int) (-min_x / mi->bin_width_x);

    mi->bin_width_y = ( max_y - min_y ) / n_ybins;
    mi->bin_offset_y = (int) (-min_y /mi->bin_width_y);

    printf("width x %f\n", mi->bin_width_x); 
    printf("width y %f\n", mi->bin_width_y); 

    printf("offset x %d\n", mi->bin_offset_x); 
    printf("offset y %d\n", mi->bin_offset_y); 

    for (int j=1;j<n_triangles;j++) {
        float xx = 0.0;
        float yy = 0.0;
        
        for (int i=0;i<3;i++) {

            int a = triangles[j*3+i];

            float x = points[a*2  ];
            float y = points[a*2+1];
//            float z = points[a*3+2];

            xx += x;
            yy += y;
        }
        xx /= 3.0;
        yy /= 3.0;

        int a = mi->bin_offset_x + (int) (xx/mi->bin_width_x) ;
        int b = mi->bin_offset_y + (int) (yy/mi->bin_width_y) ;

      //  printf("a,b,i,n %d %d %d %d \n",a,b,b * mi->n_xbins + a, n_xbins*n_ybins);
        mi->bincounts[ b * mi->n_xbins + a ]++;

    }

    for (int a = 0; a<mi->n_xbins; a++) {
        for (int b = 0; b<mi->n_ybins; b++) {
        //    printf("bincounts:  (%d,%d) %d\n",a,b, mi->bincounts[ b * mi->n_xbins + a ]);
            max_binfill = max_binfill < mi->bincounts[ b * mi->n_xbins + a ] ? mi->bincounts[ b * mi->n_xbins + a ]  : max_binfill;
        }    
    }

    mi->max_bin_triangles = max_binfill;

    mi->bins = malloc( sizeof(int)*max_binfill*mi->n_xbins*mi->n_ybins );
    memset(mi->bincounts,0,sizeof(int)*(1+n_xbins)*(1+n_ybins));



    for (int j=0;j<n_triangles;j++) {
        float xx = 0.0;
        float yy = 0.0;
        
        for (int i=0;i<3;i++) {

            int a = triangles[j*3+i];

            float x = points[a*2  ];
            float y = points[a*2+1];
//            float z = points[a*3+2];

            xx += x;
            yy += y;
        }
        xx /= 3.0;
        yy /= 3.0;

        int a = mi->bin_offset_x + (int) (xx/mi->bin_width_x) ;
        int b = mi->bin_offset_y + (int) (yy/mi->bin_width_y) ;

//        printf("a,b,i,n %d %d %d %d \n",a,b,b * mi->n_xbins + a, n_xbins*n_ybins);

        int k = mi->bincounts[ b * mi->n_xbins + a ]++;
        mi->bins[ (b * mi->n_xbins + a)*mi->max_bin_triangles + k ] = j;
//        mi->bincounts[ b * mi->n_xbins + a ]++;

    }

    return mi;

}

#endif