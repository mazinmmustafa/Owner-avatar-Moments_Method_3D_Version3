//
#include "quadl.hpp"

const real_t eps_zero=1.0E-2;

size_t max_size_t(const size_t a, const size_t b){
    return a>b ? a : b;
}

size_t min_size_t(const size_t a, const size_t b){
    return a<b ? a : b;
}

extern "C" void cgqf_f77_(int *rule, int *order, double *x, double *w);
extern "C" void dunavant_rule_wrapper_(int *rule, int *order_num, double *x, double *y, double *w);

quadl_t::quadl_t(){
}

quadl_t::~quadl_t(){
}

void quadl_t::set(const size_t N, const size_t k_max, const real_t tol){
    assert_error(this->is_allocated==false, "quad is already set");
    assert_error(tol>0.0, "invalid tolerance");
    this->N = N;
    this->k_max = k_max;
    this->x = (real_t*)calloc(this->N, sizeof(real_t));
    assert(this->x!=null);
    this->w = (real_t*)calloc(this->N, sizeof(real_t));
    assert(this->w!=null);
    int rule=1;
    int N_=this->N;
    double *x_, *w_;
    x_ = (double*)calloc(this->N, sizeof(double));
    assert(this->x!=null);
    w_ = (double*)calloc(this->N, sizeof(double));
    assert(this->w!=null);
    cgqf_f77_(&rule, &N_, x_, w_);
    this->x = (real_t*)(x_);
    this->w = (real_t*)(w_);
    this->is_allocated = true;
}

void quadl_t::unset(){
    if (this->is_allocated){
        free(x);
        free(w);
        this->N = 0;
        this->k_max = 0;
        this->is_allocated = false;
    }
}

void quadl_t::disp(){
    assert_error(this->is_allocated, "quadrature rule is unset");
    for (size_t i=0; i<this->N; i++){
        print("%2zu: %21.14E, %21.14E\n", i, this->x[i], this->w[i]);
    }
}

// 1d

complex_t quadl_t::quadl_1d(complex_t (*func)(const complex_t, void*), 
    void *args, const real_t a, const real_t b){
    real_t h_p=(b+a)/2.0;
    real_t h_m=(b-a)/2.0;
    complex_t sum=0.0;
    real_t w_n;
    real_t x_n;
    for (size_t n=0; n<this->N; n++){
        w_n = this->w[n];
        x_n = this->x[n];
        sum+=h_m*w_n*func(h_m*x_n+h_p, args);
    }
    return sum;
}

complex_t quadl_t::quadl_1d_(complex_t (*func)(const complex_t, void*), 
    void *args, const real_t a, const real_t b, size_t &k, const complex_t I_p){
    real_t m=(a+b)/2.0;
    complex_t I1=quadl_t::quadl_1d(func, args, a, m);
    complex_t I2=quadl_t::quadl_1d(func, args, m, b);
    complex_t I_n=I1+I2;
    real_t error=abs(I_n-I_p);
    if (error>this->tol*abs(I_n)&&k<this->k_max&&error>0.0){
        k++;
        I1 = quadl_t::quadl_1d_(func, args, a, m, k, I1);
        if (k>this->k_max){return I_n;}
        I2 = quadl_t::quadl_1d_(func, args, m, b, k, I2);
        if (k>this->k_max){return I_n;}
        I_n = I1+I2;
    }
    return I_n;
}

complex_t quadl_t::integral_1d(complex_t (*func)(const complex_t, void*), 
    void *args, const real_t a, const real_t b, int &flag){
    assert_error(this->is_allocated, "quadrature rule is unset");
    flag = false;
    size_t k=0;
    complex_t ans=quadl_t::quadl_1d_(func, args, a, b, k, 0.0);
    if (k>=this->k_max){flag = true;}
    return ans;
}

// 2d

complex_t quadl_t::quadl_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y){
    real_t h_x_p=(b_x+a_x)/2.0;
    real_t h_x_m=(b_x-a_x)/2.0;
    real_t h_y_p=(b_y+a_y)/2.0;
    real_t h_y_m=(b_y-a_y)/2.0;
    complex_t sum=0.0;
    real_t w_m;
    real_t x_m;
    real_t w_n;
    real_t x_n;
    for (size_t m=0; m<this->N; m++){
        w_m = this->w[m];
        x_m = this->x[m];
        for (size_t n=0; n<this->N; n++){
            w_n = this->w[n];
            x_n = this->x[n];
            sum+=h_x_m*h_y_m*w_m*w_n*func(h_x_m*x_m+h_x_p, h_y_m*x_n+h_y_p, args);
        }
    }
    return sum;
}

complex_t quadl_t::quadl_2d_(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, 
    size_t &k, const complex_t I_p){
    real_t m_x=(a_x+b_x)/2.0;
    real_t m_y=(a_y+b_y)/2.0;
    complex_t I1=quadl_t::quadl_2d(func, args, a_x, m_x, a_y, m_y);
    complex_t I2=quadl_t::quadl_2d(func, args, a_x, m_x, m_y, b_y);
    complex_t I3=quadl_t::quadl_2d(func, args, m_x, b_x, a_y, m_y);
    complex_t I4=quadl_t::quadl_2d(func, args, m_x, b_x, m_y, b_y);
    complex_t I_n=I1+I2+I3+I4;
    real_t error=abs(I_n-I_p);
    if (error>this->tol*abs(I_n)&&k<this->k_max&&error>0.0){
        k++;
        I1 = quadl_t::quadl_2d_(func, args, a_x, m_x, a_y, m_y, k, I1);
        if (k>this->k_max){return I_n;}
        I2 = quadl_t::quadl_2d_(func, args, a_x, m_x, m_y, b_y, k, I2);
        if (k>this->k_max){return I_n;}
        I3 = quadl_t::quadl_2d_(func, args, m_x, b_x, a_y, m_y, k, I3);
        if (k>this->k_max){return I_n;}
        I4 = quadl_t::quadl_2d_(func, args, m_x, b_x, m_y, b_y, k, I4);
        if (k>this->k_max){return I_n;}
        I_n = I1+I2+I3+I4; 
    }
    return I_n;
}

complex_t quadl_t::integral_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
    int &flag){
    assert_error(this->is_allocated, "quadrature rule is unset");
    flag = false;
    size_t k=0;
    complex_t ans=quadl_t::quadl_2d_(func, args, a_x, b_x, a_y, b_y, k, 0.0);
    if (k>=this->k_max){flag = true;}
    return ans;
}

// 3d

complex_t quadl_t::quadl_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
    const real_t a_z, const real_t b_z){
    real_t h_x_p=(b_x+a_x)/2.0;
    real_t h_x_m=(b_x-a_x)/2.0;
    real_t h_y_p=(b_y+a_y)/2.0;
    real_t h_y_m=(b_y-a_y)/2.0;
    real_t h_z_p=(b_z+a_z)/2.0;
    real_t h_z_m=(b_z-a_z)/2.0;
    complex_t sum=0.0;
    real_t w_p;
    real_t x_p;
    real_t w_m;
    real_t x_m;
    real_t w_n;
    real_t x_n;
    for (size_t p=0; p<this->N; p++){
        w_p = this->w[p];
        x_p = this->x[p];
        for (size_t m=0; m<this->N; m++){
            w_m = this->w[m];
            x_m = this->x[m];
            for (size_t n=0; n<this->N; n++){
                w_n = this->w[n];
                x_n = this->x[n];
                sum+=h_x_m*h_y_m*h_z_m*w_p*w_m*w_n*func(
                    h_x_m*x_m+h_x_p, h_y_m*x_n+h_y_p, h_z_p*x_p+h_z_p, args);
            }
        }
    }
    return sum;
}

complex_t quadl_t::quadl_3d_(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, 
    const real_t a_z, const real_t b_z, size_t &k, const complex_t I_p){
    real_t m_x=(a_x+b_x)/2.0;
    real_t m_y=(a_y+b_y)/2.0;
    real_t m_z=(a_z+b_z)/2.0;
    complex_t I1=quadl_t::quadl_3d(func, args, a_x, m_x, a_y, m_y, a_z, m_z);
    complex_t I2=quadl_t::quadl_3d(func, args, a_x, m_x, a_y, m_y, m_z, b_z);
    complex_t I3=quadl_t::quadl_3d(func, args, a_x, m_x, m_y, b_y, a_z, m_z);
    complex_t I4=quadl_t::quadl_3d(func, args, a_x, m_x, m_y, b_y, m_z, b_z);
    complex_t I5=quadl_t::quadl_3d(func, args, m_x, b_x, a_y, m_y, a_z, m_z);
    complex_t I6=quadl_t::quadl_3d(func, args, m_x, b_x, a_y, m_y, m_z, b_z);
    complex_t I7=quadl_t::quadl_3d(func, args, m_x, b_x, m_y, b_y, a_z, m_z);
    complex_t I8=quadl_t::quadl_3d(func, args, m_x, b_x, m_y, b_y, m_z, b_z);
    complex_t I_n=I1+I2+I3+I4+I5+I6+I7+I8;
    real_t error=abs(I_n-I_p);
    if (error>this->tol*abs(I_n)&&k<this->k_max&&error>0.0){
        k++;
        I1 = quadl_t::quadl_3d_(func, args, a_x, m_x, a_y, m_y, a_z, m_z, k, I1);
        if (k>this->k_max){return I_n;}
        I2 = quadl_t::quadl_3d_(func, args, a_x, m_x, a_y, m_y, m_z, b_z, k, I2);
        if (k>this->k_max){return I_n;}
        I3 = quadl_t::quadl_3d_(func, args, a_x, m_x, m_y, b_y, a_z, m_z, k, I3);
        if (k>this->k_max){return I_n;}
        I4 = quadl_t::quadl_3d_(func, args, a_x, m_x, m_y, b_y, m_z, b_z, k, I4);
        if (k>this->k_max){return I_n;}
        I5 = quadl_t::quadl_3d_(func, args, m_x, b_x, a_y, m_y, a_z, m_z, k, I5);
        if (k>this->k_max){return I_n;}
        I6 = quadl_t::quadl_3d_(func, args, m_x, b_x, a_y, m_y, m_z, b_z, k, I6);
        if (k>this->k_max){return I_n;}
        if (k>this->k_max){return I_n;}
        I7 = quadl_t::quadl_3d_(func, args, m_x, b_x, m_y, b_y, a_z, m_z, k, I7);
        if (k>this->k_max){return I_n;}
        I8 = quadl_t::quadl_3d_(func, args, m_x, b_x, m_y, b_y, m_z, b_z, k, I8);
        if (k>this->k_max){return I_n;}
        I_n = I1+I2+I3+I4+I5+I6+I7+I8; 
    }
    return I_n;
}

complex_t quadl_t::integral_3d(complex_t (*func)(const complex_t, const complex_t,const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
    const real_t a_z, const real_t b_z, int &flag){
    assert_error(this->is_allocated, "quadrature rule is unset");
    flag = false;
    size_t k=0;
    complex_t ans=quadl_t::quadl_3d_(func, args, a_x, b_x, a_y, b_y, a_z, b_z, k, 0.0);
    if (k>=this->k_max){flag = true;}
    return ans;
}

// 1d line domain
void quadl_domain_t::set_1d(const size_t k_max, const real_t tol){
    assert_error(tol>0.0, "invalid value for tol");
    assert_error(k_max>0, "invalid value for k_max");
    this->k_max_1d = k_max;
    this->tol_1d = tol;
}

complex_t quadl_domain_t::quadl_1d(complex_t (*func)(const complex_t, void*), 
    void *args, line_domain_t line){
    real_t L=line.length();
    complex_t sum=0.0;
    real_t alpha;
    for (size_t n=0; n<this->N_1d; n++){
        alpha = line.v1.x+(line.v2.x-line.v1.x)*this->x_1d[n];
        sum+=this->w_1d[n]*func(alpha, args);
    }
    return L*sum;
}

complex_t quadl_domain_t::quadl_1d_(complex_t (*func)(const complex_t, void*), 
    void *args, const line_domain_t line, size_t &k, const complex_t I_p){
    if (k>this->k_max_1d){return I_p;}
    vector_t<real_t> m=0.5*(line.v1+line.v2);
    line_domain_t line_1={line.v1, m};
    line_domain_t line_2={m, line.v2};
    complex_t I1=quadl_domain_t::quadl_1d(func, args, line_1);
    complex_t I2=quadl_domain_t::quadl_1d(func, args, line_2);
    complex_t I_n=I1+I2;
    real_t error=abs(I_n-I_p);
    size_t k1=k, k2=k;
    if (error>this->tol_1d*abs(I_n)&&error>0.0&&abs(I_n)>eps_zero){
        k++;
        I1 = quadl_domain_t::quadl_1d_(func, args, line_1, ++k1, I1);
        I2 = quadl_domain_t::quadl_1d_(func, args, line_2, ++k2, I2);
        I_n = I1+I2;
        k = max_size_t(k1, k2);
    }
    return I_n;
}

complex_t quadl_domain_t::integral_1d(complex_t (*func)(const complex_t, void*), 
    void *args, const line_domain_t line, int &flag){
    flag = false;
    size_t k=0;
    complex_t ans=quadl_domain_t::quadl_1d_(func, args, line, k, 0.0);
    if (k>=this->k_max_1d){flag = true;}
    return ans;
}

// 2d domain triangle

void quadl_domain_t::set_2d(const size_t k_max, const real_t tol){
    assert_error(tol>0.0, "invalid value for tol");
    assert_error(k_max>0, "invalid value for k_max");
    this->k_max_2d = k_max;
    this->tol_2d = tol;
}

complex_t quadl_domain_t::quadl_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, triangle_domain_t triangle){
    real_t A=triangle.area();
    complex_t sum=0.0;
    real_t alpha, beta;
    for (size_t n=0; n<this->N_2d; n++){
        alpha = triangle.v1.x+
                (triangle.v2.x-triangle.v1.x)*this->x_2d[n]+
                (triangle.v3.x-triangle.v1.x)*this->y_2d[n];
        beta  = triangle.v1.y+
                (triangle.v2.y-triangle.v1.y)*this->x_2d[n]+
                (triangle.v3.y-triangle.v1.y)*this->y_2d[n];
        sum+=this->w_2d[n]*func(alpha, beta, args);
    }
    return A*sum;
}

complex_t quadl_domain_t::quadl_2d_(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, const triangle_domain_t triangle, size_t &k, const complex_t I_p){
    if (k>this->k_max_2d){return I_p;}
    //
    triangle_domain_t triangle_1, triangle_2;
    vector_t<real_t> m23=0.5*(triangle.v2+triangle.v3);
    //
    triangle_1.v1 = triangle.v1; 
    triangle_1.v2 = triangle.v2; 
    triangle_1.v3 = m23;
    //
    triangle_2.v1 = triangle.v1; 
    triangle_2.v2 = m23;
    triangle_2.v3 = triangle.v3; 
    //
    complex_t I1=quadl_domain_t::quadl_2d(func, args, triangle_1);
    complex_t I2=quadl_domain_t::quadl_2d(func, args, triangle_2);
    complex_t I_n=I1+I2;
    real_t error=abs(I_n-I_p);
    size_t k1=k, k2=k;
    if (error>this->tol_1d*abs(I_n)&&error>0.0&&abs(I_n)>eps_zero){
        I1 = quadl_domain_t::quadl_2d_(func, args, triangle_1, ++k1, I1);
        I2 = quadl_domain_t::quadl_2d_(func, args, triangle_2, ++k2, I2);
        I_n = I1+I2;
        k = max_size_t(k1, k2);
    }
    return I_n;
}

complex_t quadl_domain_t::integral_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, const triangle_domain_t triangle, int &flag){
    flag = false;
    size_t k=0;
    complex_t ans=quadl_domain_t::quadl_2d_(func, args, triangle, k, 0.0);
    if (k>=this->k_max_2d){flag = true; print(ans);}
    return ans;
}

// 3d domain tetrahedron

void quadl_domain_t::set_3d(const size_t k_max, const real_t tol){
    assert_error(tol>0.0, "invalid value for tol");
    assert_error(k_max>0, "invalid value for k_max");
    this->k_max_3d = k_max;
    this->tol_3d = tol;
}

complex_t quadl_domain_t::quadl_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
    void *args, tetrahedron_domain_t tetrahedron){
    real_t V=tetrahedron.volume();
    complex_t sum=0.0;
    real_t alpha, beta, gamma;
    for (size_t n=0; n<this->N_3d; n++){
        alpha =  tetrahedron.v1.x+
                (tetrahedron.v2.x-tetrahedron.v1.x)*x_3d[n]+
                (tetrahedron.v3.x-tetrahedron.v1.x)*y_3d[n]+
                (tetrahedron.v4.x-tetrahedron.v1.x)*z_3d[n];
        beta  =  tetrahedron.v1.y+
                (tetrahedron.v2.y-tetrahedron.v1.y)*x_3d[n]+
                (tetrahedron.v3.y-tetrahedron.v1.y)*y_3d[n]+
                (tetrahedron.v4.y-tetrahedron.v1.y)*z_3d[n];
        gamma =  tetrahedron.v1.z+
                (tetrahedron.v2.z-tetrahedron.v1.z)*x_3d[n]+
                (tetrahedron.v3.z-tetrahedron.v1.z)*y_3d[n]+
                (tetrahedron.v4.z-tetrahedron.v1.z)*z_3d[n];
        sum+=w_3d[n]*func(alpha, beta, gamma, args);
    }
    return V*sum;
}

complex_t quadl_domain_t::quadl_3d_(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
    void *args, const tetrahedron_domain_t tetrahedron, size_t &k, const complex_t I_p){
    if (k>this->k_max_3d){return I_p;}
    tetrahedron_domain_t tetrahedron_1, tetrahedron_2, tetrahedron_3, tetrahedron_4;
    vector_t<real_t> m23=0.5*(tetrahedron.v2+tetrahedron.v3);
    vector_t<real_t> m24=0.5*(tetrahedron.v2+tetrahedron.v4);
    vector_t<real_t> m34=0.5*(tetrahedron.v3+tetrahedron.v4);
    //
    tetrahedron_1.v1 = tetrahedron.v1; 
    tetrahedron_1.v2 = tetrahedron.v4; 
    tetrahedron_1.v3 = m24; 
    tetrahedron_1.v4 = m34; 
    //
    tetrahedron_2.v1 = tetrahedron.v1; 
    tetrahedron_2.v2 = tetrahedron.v2; 
    tetrahedron_2.v3 = m23; 
    tetrahedron_2.v4 = m24; 
    //
    tetrahedron_3.v1 = tetrahedron.v1; 
    tetrahedron_3.v2 = tetrahedron.v3; 
    tetrahedron_3.v3 = m34; 
    tetrahedron_3.v4 = m23; 
    //
    tetrahedron_4.v1 = tetrahedron.v1; 
    tetrahedron_4.v2 = m24; 
    tetrahedron_4.v3 = m23;  
    tetrahedron_4.v4 = m34; 
    //
    complex_t I1=quadl_domain_t::quadl_3d(func, args, tetrahedron_1);
    complex_t I2=quadl_domain_t::quadl_3d(func, args, tetrahedron_2);
    complex_t I3=quadl_domain_t::quadl_3d(func, args, tetrahedron_3);
    complex_t I4=quadl_domain_t::quadl_3d(func, args, tetrahedron_4);
    complex_t I_n=I1+I2+I3+I4;
    real_t error=abs(I_n-I_p);
    size_t k1=k, k2=k, k3=k, k4=k;
    if (error>this->tol_3d*abs(I_n)&&error>0.0&&abs(I_n)>eps_zero){
        k++;
        I1 = quadl_domain_t::quadl_3d_(func, args, tetrahedron_1, ++k1, I1);
        I2 = quadl_domain_t::quadl_3d_(func, args, tetrahedron_2, ++k2, I2);
        I3 = quadl_domain_t::quadl_3d_(func, args, tetrahedron_3, ++k3, I3);
        I4 = quadl_domain_t::quadl_3d_(func, args, tetrahedron_4, ++k4, I4);
        I_n = I1+I2+I3+I4;
        k = max_size_t(max_size_t(max_size_t(k1, k2), k3), k4);
    }
    return I_n;
}

complex_t quadl_domain_t::integral_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
    void *args, const tetrahedron_domain_t tetrahedron, int &flag){
    flag = false;
    size_t k=0;
    complex_t ans=quadl_domain_t::quadl_3d_(func, args, tetrahedron, k, 0.0);
    if (k>=this->k_max_3d){flag = true;}
    return ans;
}
