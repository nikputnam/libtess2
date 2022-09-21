

void normalize(float* n, float* r);

void cross_product(float* a, float* b, float* n);

void subtract(float* r, float* a, float* b);

void add(float* r, float* a, float* b);

void subtract2(float* r, float* a, float* b);

void add2(float* r, float* a, float* b);
void set2(float* r, float* a) ;

void set(float* r, float* a) ;

void weighted_sum2(float* r, float w1, float w2, float* v1, float* v2) ;	
void weighted_sum2_4(float* r, float w1, float w2, float w3, float w4, float* v1, float* v2, float* v3, float* v4) ;	

float dist22(float* a, float* b) ;
float dist2(float* a, float* b) ;

float norm22(float* n);
float norm(float* n);

int vequal2(float* a, float* b);
int vequal3(float* a, float* b);

void scale(float* p, float a) ;

float dot2(float* x, float* y);

void diff(float* x, float* y, float* z) ;

float dot(float* x, float* y) ;
