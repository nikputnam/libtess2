




void catmullrom2(float t, float* points, int npoints, float* x);
void catmullrom2_tangent(float t, float* points, int npoints, float* x);
void catmullrom2_segment(float t, float* p0, float* p1, float* p2, float* p3, float alpha, float* x) ;
void catmullrom2_segment_tangent(float t, float* p0, float* p1, float* p2, float* p3, float alpha, float* x) ;
