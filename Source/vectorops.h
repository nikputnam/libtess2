

void normalize(float* n, float* r);

void cross_product(float* a, float* b, float* n);

void subtract(float* r, float* a, float* b);

void add(float* r, float* a, float* b);

void subtract2(float* r, float* a, float* b);

void add2(float* r, float* a, float* b);
void set2(float* r, float* a) ;

void weighted_sum2(float* r, float w1, float w2, float* v1, float* v2) ;	

float dist22(float* a, float* b) ;
float dist2(float* a, float* b) ;

float norm22(float* n);
float norm(float* n);

int vequal2(float* a, float* b);
int vequal3(float* a, float* b);