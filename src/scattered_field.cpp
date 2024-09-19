//
#include "scattered_field.hpp"

// 1d

void integrand_L1_1d(basis_1d_t b_m, const vector_t<real_t> p, 
    real_t &I_m, real_t &I_p, const real_t a){
    projection_1d_para para;
    real_t ans;
    // m
    para = prjection_1d(b_m.r_m, b_m.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    I_m = para.P_0 == 0.0 ? 0.0 : ans;
    // p
    para = prjection_1d(b_m.e[0], b_m.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    I_p = para.P_0 == 0.0 ? 0.0 : ans;
}

void integrand_L2_1d(basis_1d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p, const real_t a){
    projection_1d_para para;
    vector_t<real_t> ans;
    // m
    para = prjection_1d(b_m.r_m, b_m.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_m = para.P_0 == 0.0 ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
    // p
    para = prjection_1d(b_m.e[0], b_m.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (para.P_p-para.P_m)*para.l_unit;
    I_p = para.P_0 == 0.0 ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
}

void integrand_L3_1d(basis_1d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p, const real_t a){
    projection_1d_para para;
    vector_t<real_t> ans;
    // m
    para = prjection_1d(b_m.r_m, b_m.e[0], p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
         -(para.l_p/para.P_p-para.l_m/para.P_m)*para.P_0_unit/para.P_0;
    I_m = para.P_0 == 0.0 ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
    // p
    para = prjection_1d(b_m.e[0], b_m.r_p, p);
    para.P_p = sqrt(para.P_p*para.P_p+a*a);
    para.P_m = sqrt(para.P_m*para.P_m+a*a);
    ans = (1.0/para.P_p-1.0/para.P_m)*para.l_unit
         -(para.l_p/para.P_p-para.l_m/para.P_m)*para.P_0_unit/para.P_0;
    I_p = para.P_0 == 0.0 ? vector_t<real_t>(0.0, 0.0, 0.0) : ans;
}

complex_t E_1d_singular_integrand_1(const complex_t alpha, void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    const complex_t k=args->k;
    const real_t a=args->a;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    R_m = sqrt(R_m*R_m+a*a);
    R_p = sqrt(R_p*R_p+a*a);
    I_m = -j*k*exp(-j*k*R_m/2.0)*sinc(k*R_m/2.0)*(+1.0*alpha*b_m.L_m[0]*args->unit_vector);
    I_p = -j*k*exp(-j*k*R_p/2.0)*sinc(k*R_p/2.0)*(-1.0*alpha*b_m.L_p[0]*args->unit_vector);
    return (I_m+I_p)/(4.0*pi);
}

complex_t E_1d_singular_integrand_2(const complex_t alpha, void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    const complex_t k=args->k;
    const real_t a=args->a;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    R_m = sqrt(R_m*R_m+a*a);
    R_p = sqrt(R_p*R_p+a*a);
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*(unit(R_m_vector-r)*args->unit_vector);
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*(unit(R_p_vector-r)*args->unit_vector);
    return (I_m-I_p)/(4.0*pi);
}

complex_t E_1d_integral_1(void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    real_t a=args->a;
    const vector_t<real_t> r=args->r;
    real_t I_m, I_p;
    integrand_L1_1d(b_m, r, I_m, I_p, a);
    projection_1d_para para_m = prjection_1d(b_m.r_m, b_m.e[0], r);
    projection_1d_para para_p = prjection_1d(b_m.e[0], b_m.r_p, r);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*(para_m.l_m*para_m.l_unit)*I_m/mag(b_m.L_m[0]);
    ans+= +1.0*args->unit_vector*(para_p.l_p*para_p.l_unit)*I_p/mag(b_m.L_p[0]);
    return ans/(4.0*pi);
}

complex_t E_1d_integral_2(void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    real_t a=args->a;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L2_1d(b_m, r, I_m, I_p, a);
    complex_t ans=0.0;
    ans+= +1.0*args->unit_vector*I_m/mag(b_m.L_m[0]);
    ans+= -1.0*args->unit_vector*I_p/mag(b_m.L_p[0]);
    return ans/(4.0*pi);
}

complex_t E_1d_integral_3(void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    real_t a=args->a;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_1d(b_m, r, I_m, I_p, a);
    complex_t ans=0.0;
    ans+= +1.0*args->unit_vector*I_m/mag(b_m.L_m[0]);
    ans+= -1.0*args->unit_vector*I_p/mag(b_m.L_p[0]);
    return ans/(4.0*pi);

}

complex_t compute_E_1d(const basis_1d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, 
    const real_t a, quadl_domain_t quadl){
    scattered_field_args_1d_t args;
    const complex_t j=complex_t(0.0, 1.0);
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.a = a;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2, I3, I4, I5;
    int_t flag;
    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0)};
    //
    I1 = quadl.integral_1d(E_1d_singular_integrand_1, &args, line, flag);
    assert_error(!flag, "no convergence");
    I2 = E_1d_integral_1(&args);
    I3 = E_1d_integral_2(&args);
    I4 = quadl.integral_1d(E_1d_singular_integrand_2, &args, line, flag);
    assert_error(!flag, "no convergence");
    I5 = E_1d_integral_3(&args);
    return -j*k*eta*(I1+I2+I3)+j*(eta/k)*(I4+I5);
}

complex_t H_1d_singular_integrand_1(const complex_t alpha, void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    const complex_t k=args->k;
    const real_t a=args->a;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    R_m = sqrt(R_m*R_m+a*a);
    R_p = sqrt(R_p*R_p+a*a);
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*
        (alpha*(+1.0*b_m.L_m[0]^unit(R_m_vector-r))*args->unit_vector);
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*
        (alpha*(-1.0*b_m.L_p[0]^unit(R_p_vector-r))*args->unit_vector);
    return (I_m+I_p)/(4.0*pi);
}

complex_t H_1d_integral_1(void *args_){
    scattered_field_args_1d_t *args=(scattered_field_args_1d_t*)args_;
    basis_1d_t b_m=args->b_m;
    real_t a=args->a;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_1d(b_m, r, I_m, I_p, a);
    projection_1d_para para_m=prjection_1d(b_m.r_m, b_m.e[0], r);
    projection_1d_para para_p=prjection_1d(b_m.e[0], b_m.r_p, r);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*((para_m.l_m*para_m.l_unit)^I_m)/mag(b_m.L_m[0]);
    ans+= +1.0*args->unit_vector*((para_p.l_p*para_p.l_unit)^I_p)/mag(b_m.L_p[0]);
    return ans/(4.0*pi);
}

complex_t compute_H_1d(const basis_1d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, 
    const real_t a, quadl_domain_t quadl){
    scattered_field_args_1d_t args;
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.a = a;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2;
    int_t flag;
    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0)};
    //
    I1 = quadl.integral_1d(H_1d_singular_integrand_1, &args, line, flag);
    assert_error(!flag, "no convergence");
    I2 = H_1d_integral_1(&args);
    return I1+I2;
}

// 2d 

void integrand_L1_2d(basis_2d_t b_m, const vector_t<real_t> p, 
    real_t &I_m, real_t &I_p){
    projection_2d_para para;
    // m
    para = prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], p);
    I_m = 0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
        B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
        C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
        D = isnan(D) ? 0.0 : D; D = isinf(D) ? 0.0 : D;
        I_m+=A*(B-abs(para.d)*(C-D));
    }
    // p
    para = prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], p);
    I_p = 0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
        B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
        C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
        D = isnan(D) ? 0.0 : D; D = isinf(D) ? 0.0 : D;
        I_p+=A*(B-abs(para.d)*(C-D));
    }
}

void integrand_L2_2d(basis_2d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p){
    projection_2d_para para;
    // m
    para = prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], p);
    I_m = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
        B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
        I_m = I_m+0.5*(A+B)*para.u[i];
    }
    // p
    para = prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], p);
    I_p = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
        B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
        I_p = I_p+0.5*(A+B)*para.u[i];
    }
}

void integrand_L3_2d(basis_2d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p){
    projection_2d_para para;
    // m
    para = prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], p);
    I_m = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
        B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
        C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
        I_m = I_m+A*para.u[i]+sign(para.d)*(B-C)*(para.u[i]*para.para_1d[i].P_0_unit)*para.n;
    }
    // p
    para = prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], p);
    I_p = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
        B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
        C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
        I_p = I_p+A*para.u[i]+sign(para.d)*(B-C)*(para.u[i]*para.para_1d[i].P_0_unit)*para.n;
    }
}

complex_t E_2d_singular_integrand_1(const complex_t alpha, const complex_t beta, void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    I_m = -j*k*exp(-j*k*R_m/2.0)*sinc(k*R_m/2.0)*b_m.L*(+1.0*(alpha*b_m.L_m[0]+beta*b_m.L_m[1])*args->unit_vector);
    I_p = -j*k*exp(-j*k*R_p/2.0)*sinc(k*R_p/2.0)*b_m.L*(-1.0*(alpha*b_m.L_p[1]+beta*b_m.L_p[0])*args->unit_vector);
    return (I_m+I_p)/(4.0*pi);
}

complex_t E_2d_singular_integrand_2(const complex_t alpha, const complex_t beta, void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*(unit(R_m_vector-r)*args->unit_vector)*2.0*b_m.L;
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*(unit(R_p_vector-r)*args->unit_vector)*2.0*b_m.L;
    return (I_m-I_p)/(4.0*pi);
}

complex_t E_2d_integral_1(void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    real_t I_m, I_p;
    integrand_L1_2d(b_m, r, I_m, I_p);
    projection_2d_para para_m = prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], r);
    projection_2d_para para_p = prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], r);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*(b_m.r_m-para_m.p_0)*I_m*b_m.L/(2.0*b_m.A_m[0]);
    ans+= +1.0*args->unit_vector*(b_m.r_p-para_p.p_0)*I_p*b_m.L/(2.0*b_m.A_p[0]);
    return ans/(4.0*pi);
}

complex_t E_2d_integral_2(void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L2_2d(b_m, r, I_m, I_p);
    complex_t ans=0.0;
    ans+= +1.0*args->unit_vector*I_m*b_m.L/(2.0*b_m.A_m[0]);
    ans+= -1.0*args->unit_vector*I_p*b_m.L/(2.0*b_m.A_p[0]);
    return ans/(4.0*pi);
}

complex_t E_2d_integral_3(void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_2d(b_m, r, I_m, I_p);
    complex_t ans=0.0;
    ans+= +1.0*args->unit_vector*I_m*b_m.L/b_m.A_m[0];
    ans+= -1.0*args->unit_vector*I_p*b_m.L/b_m.A_p[0];
    return ans/(4.0*pi);
}

complex_t compute_E_2d(const basis_2d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, quadl_domain_t quadl){
    scattered_field_args_2d_t args;
    const complex_t j=complex_t(0.0, 1.0);
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2, I3, I4, I5;
    int_t flag;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0), vector_t<real_t>(0.0, +1.0, 0.0)};
    //
    I1 = quadl.integral_2d(E_2d_singular_integrand_1, &args, triangle, flag);
    assert_error(!flag, "no convergence");
    I2 = E_2d_integral_1(&args);
    I3 = E_2d_integral_2(&args);
    I4 = quadl.integral_2d(E_2d_singular_integrand_2, &args, triangle, flag);
    assert_error(!flag, "no convergence");
    I5 = E_2d_integral_3(&args);
    return -j*k*eta*(I1+I2+I3)+j*(eta/k)*(I4+I5);
}

complex_t H_2d_singular_integrand_1(const complex_t alpha, const complex_t beta, void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[1]+real(beta)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*
        ((+1.0*(alpha*b_m.L_m[0]+beta*b_m.L_m[1])^unit(R_m_vector-r))*args->unit_vector)*b_m.L;
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*
        ((-1.0*(alpha*b_m.L_p[1]+beta*b_m.L_p[0])^unit(R_p_vector-r))*args->unit_vector)*b_m.L;
    return (I_m+I_p)/(4.0*pi);
}

complex_t H_2d_integral_1(void *args_){
    scattered_field_args_2d_t *args=(scattered_field_args_2d_t*)args_;
    basis_2d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_2d(b_m, r, I_m, I_p);
    projection_2d_para para_m=prjection_2d(b_m.r_m, b_m.e[0], b_m.e[1], r);
    projection_2d_para para_p=prjection_2d(b_m.r_p, b_m.e[1], b_m.e[0], r);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*((b_m.r_m-para_m.p_0)^I_m)*b_m.L/(2.0*b_m.A_m[0]);
    ans+= +1.0*args->unit_vector*((b_m.r_p-para_p.p_0)^I_p)*b_m.L/(2.0*b_m.A_p[0]);
    return ans/(4.0*pi);
}

complex_t compute_H_2d(const basis_2d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, quadl_domain_t quadl){
    scattered_field_args_2d_t args;
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2;
    int_t flag;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(+1.0, 0.0, 0.0), vector_t<real_t>(0.0, +1.0, 0.0)};
    //
    I1 = quadl.integral_2d(H_2d_singular_integrand_1, &args, triangle, flag);
    assert_error(!flag, "no convergence");
    I2 = H_2d_integral_1(&args);
    return I1+I2;
}

// 3d 

void integrand_L1_3d(basis_3d_t b_m, const vector_t<real_t> p, 
    real_t &I_m, real_t &I_p){
    projection_3d_para para;
    // m
    para = prjection_3d(b_m.r_m, b_m.e[0], b_m.e[1], b_m.e[2], p);
    I_m = 0.0;
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
            B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
            C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
            D = isnan(D) ? 0.0 : D; D = isinf(D) ? 0.0 : D;
            I_m = I_m+0.5*para.para_2d[j].d*A*(abs(para.para_2d[j].d)*(C-D)-B);
        }
    }
    // p
    para = prjection_3d(b_m.r_p, b_m.e[2], b_m.e[1], b_m.e[0], p);
    I_p = 0.0;
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
            B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
            C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
            D = isnan(D) ? 0.0 : D; D = isinf(D) ? 0.0 : D;
            I_p = I_p+0.5*para.para_2d[j].d*A*(abs(para.para_2d[j].d)*(C-D)-B);
        }
    }
}

void integrand_L2_3d(basis_3d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p){
    projection_3d_para para;
    // m
    para = prjection_3d(b_m.r_m, b_m.e[0], b_m.e[1], b_m.e[2], p);
    I_m = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=0.5*para.para_2d[j].para_1d[i].P_0*(pow(para.para_2d[j].R_0[i], 2.0)+2.0*pow(para.para_2d[j].d, 2.0));
            real_t C=
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t D=0.5*para.para_2d[j].para_1d[i].P_0*(
                    para.para_2d[j].R_p[i]*para.para_2d[j].para_1d[i].l_p-
                    para.para_2d[j].R_m[i]*para.para_2d[j].para_1d[i].l_m
            );
            real_t E=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t F=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
            B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
            C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
            D = isnan(D) ? 0.0 : D; D = isinf(D) ? 0.0 : D;
            E = isnan(E) ? 0.0 : E; E = isinf(E) ? 0.0 : E;
            I_m = I_m+(1.0/3.0)*A*(B*C+D-pow(abs(para.para_2d[j].d), 3.0)*(E-F))*para.para_2d[j].n;
        }
    }
    // p
    para = prjection_3d(b_m.r_p, b_m.e[2], b_m.e[1], b_m.e[0], p);
    I_p = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=0.5*para.para_2d[j].para_1d[i].P_0*(pow(para.para_2d[j].R_0[i], 2.0)+2.0*pow(para.para_2d[j].d, 2.0));
            real_t C=
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t D=0.5*para.para_2d[j].para_1d[i].P_0*(
                    para.para_2d[j].R_p[i]*para.para_2d[j].para_1d[i].l_p-
                    para.para_2d[j].R_m[i]*para.para_2d[j].para_1d[i].l_m
            );
            real_t E=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t F=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
            B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
            C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
            D = isnan(D) ? 0.0 : D; D = isinf(D) ? 0.0 : D;
            E = isnan(E) ? 0.0 : E; E = isinf(E) ? 0.0 : E;
            I_p = I_p+(1.0/3.0)*A*(B*C+D-pow(abs(para.para_2d[j].d), 3.0)*(E-F))*para.para_2d[j].n;
        }
    }
}

void integrand_L3_3d(basis_3d_t b_m, const vector_t<real_t> p, 
    vector_t<real_t> &I_m, vector_t<real_t> &I_p){
    projection_3d_para para;
    // m
    para = prjection_3d(b_m.r_m, b_m.e[0], b_m.e[1], b_m.e[2], p);
    I_m = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
            B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
            C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
            D = isnan(D) ? 0.0 : D; D = isinf(D) ? 0.0 : D;
            I_m = I_m+A*(B-abs(para.para_2d[j].d)*(C-D))*para.para_2d[j].n;
        }
    }
    // p
    para = prjection_3d(b_m.r_p, b_m.e[2], b_m.e[1], b_m.e[0], p);
    I_p = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t j=0; j<4; j++){
        for (size_t i=0; i<3; i++){
            real_t A=para.para_2d[j].para_1d[i].P_0_unit*para.para_2d[j].u[i];
            real_t B=para.para_2d[j].para_1d[i].P_0*
                log((para.para_2d[j].R_p[i]+para.para_2d[j].para_1d[i].l_p)/
                    (para.para_2d[j].R_m[i]+para.para_2d[j].para_1d[i].l_m));
            real_t C=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_p, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_p[i]);
            real_t D=atan2(para.para_2d[j].para_1d[i].P_0*para.para_2d[j].para_1d[i].l_m, 
                            pow(para.para_2d[j].R_0[i], 2.0)+abs(para.para_2d[j].d)*para.para_2d[j].R_m[i]);
            A = isnan(A) ? 0.0 : A; A = isinf(A) ? 0.0 : A;
            B = isnan(B) ? 0.0 : B; B = isinf(B) ? 0.0 : B;
            C = isnan(C) ? 0.0 : C; C = isinf(C) ? 0.0 : C;
            D = isnan(D) ? 0.0 : D; D = isinf(D) ? 0.0 : D;
            I_p = I_p+A*(B-abs(para.para_2d[j].d)*(C-D))*para.para_2d[j].n;
        }
    }
}

complex_t E_3d_singular_integrand_1(const complex_t alpha, const complex_t beta, const complex_t gamma, void *args_){
    scattered_field_args_3d_t *args=(scattered_field_args_3d_t*)args_;
    basis_3d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]+real(gamma)*b_m.L_m[2];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[2]+real(beta)*b_m.L_p[1]+real(gamma)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    complex_t kappa_m=(b_m.eps_m-1.0)/b_m.eps_m;
    complex_t kappa_p=(b_m.eps_p-1.0)/b_m.eps_p;
    I_m = -j*k*exp(-j*k*R_m/2.0)*sinc(k*R_m/2.0)*b_m.A*(+1.0*(alpha*b_m.L_m[0]+beta*b_m.L_m[1]+gamma*b_m.L_m[2])*args->unit_vector);
    I_p = -j*k*exp(-j*k*R_p/2.0)*sinc(k*R_p/2.0)*b_m.A*(-1.0*(alpha*b_m.L_p[2]+beta*b_m.L_p[1]+gamma*b_m.L_p[0])*args->unit_vector);
    return 2.0*(I_m*kappa_m+I_p*kappa_p)/(4.0*pi);
}

complex_t E_3d_singular_integrand_2(const complex_t alpha, const complex_t beta, const complex_t gamma, void *args_){
    scattered_field_args_3d_t *args=(scattered_field_args_3d_t*)args_;
    basis_3d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]+real(gamma)*b_m.L_m[2];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[2]+real(beta)*b_m.L_p[1]+real(gamma)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    complex_t kappa_m=(b_m.eps_m-1.0)/b_m.eps_m;
    complex_t kappa_p=(b_m.eps_p-1.0)/b_m.eps_p;
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*(unit(R_m_vector-r)*args->unit_vector)*6.0*b_m.A;
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*(unit(R_p_vector-r)*args->unit_vector)*6.0*b_m.A;
    return (I_m*kappa_m-I_p*kappa_p)/(4.0*pi);
}

complex_t E_3d_integral_1(void *args_){
    scattered_field_args_3d_t *args=(scattered_field_args_3d_t*)args_;
    basis_3d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    real_t I_m, I_p;
    integrand_L1_3d(b_m, r, I_m, I_p);
    complex_t ans=0.0;
    complex_t kappa_m=(b_m.eps_m-1.0)/b_m.eps_m;
    complex_t kappa_p=(b_m.eps_p-1.0)/b_m.eps_p;
    ans+= -1.0*args->unit_vector*(b_m.r_m-r)*I_m*kappa_m*b_m.A/(3.0*b_m.V_m);
    ans+= +1.0*args->unit_vector*(b_m.r_p-r)*I_p*kappa_p*b_m.A/(3.0*b_m.V_p);
    return ans/(4.0*pi);
}

complex_t E_3d_integral_2(void *args_){
    scattered_field_args_3d_t *args=(scattered_field_args_3d_t*)args_;
    basis_3d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L2_3d(b_m, r, I_m, I_p);
    complex_t ans=0.0;
    complex_t kappa_m=(b_m.eps_m-1.0)/b_m.eps_m;
    complex_t kappa_p=(b_m.eps_p-1.0)/b_m.eps_p;
    ans+= +1.0*args->unit_vector*I_m*kappa_m*b_m.A/(3.0*b_m.V_m);
    ans+= -1.0*args->unit_vector*I_p*kappa_p*b_m.A/(3.0*b_m.V_p);
    return ans/(4.0*pi);
}

complex_t E_3d_integral_3(void *args_){
    scattered_field_args_3d_t *args=(scattered_field_args_3d_t*)args_;
    basis_3d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_3d(b_m, r, I_m, I_p);
    complex_t ans=0.0;
    complex_t kappa_m=(b_m.eps_m-1.0)/b_m.eps_m;
    complex_t kappa_p=(b_m.eps_p-1.0)/b_m.eps_p;
    ans+= +1.0*args->unit_vector*I_m*kappa_m*b_m.A/b_m.V_m;
    ans+= -1.0*args->unit_vector*I_p*kappa_p*b_m.A/b_m.V_p;
    return ans/(4.0*pi);
}

complex_t E_3d_integral_4(void *args_){
    scattered_field_args_3d_t *args=(scattered_field_args_3d_t*)args_;
    basis_3d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> rho_m=b_m.r_m+b_m.L_m[0]/3.0+b_m.L_m[1]/3.0+b_m.L_m[2]/3.0;
    vector_t<real_t> R_vector;
    R_vector = r-rho_m;
    real_t R=mag(R_vector);
    complex_t ans=0.0;  
    real_t factor=b_m.nA*(b_m.L_m[0]/3.0+b_m.L_m[1]/3.0+b_m.L_m[2]/3.0);
    const complex_t k=args->k;
    const complex_t j=complex_t(0.0, 1.0);
    complex_t I;
    I = -j*k*(1.0+1.0/(j*k*R))*exp(-j*k*R)/R;
    ans+= +1.0*factor*(I*unit(R_vector)*args->unit_vector)/(3.0*b_m.V_m);
    return 2.0*b_m.A*ans/(4.0*pi);
}

complex_t compute_E_3d(const basis_3d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, quadl_domain_t quadl){
    scattered_field_args_3d_t args;
    const complex_t j=complex_t(0.0, 1.0);
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2, I3, I4, I5;
    int_t flag;
    tetrahedron_domain_t tetrahedron={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), 
            vector_t<real_t>(0.0, 1.0, 0.0), vector_t<real_t>(0.0, 0.0, 1.0)};
    //
    I1 = quadl.integral_3d(E_3d_singular_integrand_1, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    I2 = E_3d_integral_1(&args);
    I3 = E_3d_integral_2(&args);
    I4 = quadl.integral_3d(E_3d_singular_integrand_2, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    I5 = E_3d_integral_3(&args);
    //
    complex_t I6=0.0;
    I6 = E_3d_integral_4(&args);
    complex_t kappa_m=(b_m.eps_m-1.0)/b_m.eps_m;
    complex_t kappa_p=(b_m.eps_p-1.0)/b_m.eps_p;
    return -j*k*eta*(I1+I2+I3)+j*(eta/k)*(I4+I5)+j*(eta/k)*(kappa_p-kappa_m)*I6;
}

complex_t H_3d_singular_integrand_1(const complex_t alpha, const complex_t beta, const complex_t gamma, void *args_){
    scattered_field_args_3d_t *args=(scattered_field_args_3d_t*)args_;
    basis_3d_t b_m=args->b_m;
    const complex_t k=args->k;
    const vector_t<real_t> r=args->r;
    complex_t I_m, I_p;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_m_vector=b_m.r_m+real(alpha)*b_m.L_m[0]+real(beta)*b_m.L_m[1]+real(gamma)*b_m.L_m[2];
    vector_t<real_t> R_p_vector=b_m.r_p+real(alpha)*b_m.L_p[2]+real(beta)*b_m.L_p[1]+real(gamma)*b_m.L_p[0];
    real_t R_m, R_p;
    R_m = mag(R_m_vector-r);
    R_p = mag(R_p_vector-r);
    I_m = -0.5*k*k*exp(-j*k*R_m/2.0)*(sinc(k*R_m/2.0)+j*sinc_dx(k*R_m/2.0))*
        ((+1.0*(alpha*b_m.L_m[0]+beta*b_m.L_m[1]+gamma*b_m.L_m[2])^unit(R_m_vector-r))*args->unit_vector)*b_m.A;
    I_p = -0.5*k*k*exp(-j*k*R_p/2.0)*(sinc(k*R_p/2.0)+j*sinc_dx(k*R_p/2.0))*
        ((-1.0*(alpha*b_m.L_p[2]+beta*b_m.L_p[1]+gamma*b_m.L_p[0])^unit(R_p_vector-r))*args->unit_vector)*b_m.A;
    return 2.0*(I_m+I_p)/(4.0*pi);
}

complex_t H_3d_integral_1(void *args_){
    scattered_field_args_3d_t *args=(scattered_field_args_3d_t*)args_;
    basis_3d_t b_m=args->b_m;
    const vector_t<real_t> r=args->r;
    vector_t<real_t> I_m, I_p;
    integrand_L3_3d(b_m, r, I_m, I_p);
    complex_t ans=0.0;
    ans+= -1.0*args->unit_vector*((b_m.r_m-r)^I_m)*b_m.A/(3.0*b_m.V_m);
    ans+= +1.0*args->unit_vector*((b_m.r_p-r)^I_p)*b_m.A/(3.0*b_m.V_p);
    return ans/(4.0*pi);
}

complex_t compute_H_3d(const basis_3d_t b_m, const vector_t<real_t> r, const vector_t<real_t> unit_vector, 
    const complex_t k, const complex_t eta, quadl_domain_t quadl){
    scattered_field_args_3d_t args;
    args.r = r;
    args.unit_vector = unit_vector;
    args.b_m = b_m;
    args.k = k;
    args.eta = eta;
    complex_t I1, I2;
    int_t flag;
    tetrahedron_domain_t tetrahedron={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), 
            vector_t<real_t>(0.0, 1.0, 0.0), vector_t<real_t>(0.0, 0.0, 1.0)};
    //
    I1 = quadl.integral_3d(H_3d_singular_integrand_1, &args, tetrahedron, flag);
    assert_error(!flag, "no convergence");
    I2 = H_3d_integral_1(&args);
    return I1+I2;
}