//
#include "testbench_2d.hpp"


void test_engine_2d_2d(){

    const complex_t eta=sqrt(mu_0/eps_0);

    stopwatch_t T;
    quadl_domain_t quadl;   
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set_2d(k_max, tol);
    complex_t ans, psi, phi;
    int_t flag;
    vector_t<real_t> r_m_m, r_m_p, r_n_m, r_n_p, e_m_1, e_m_2, e_n_1, e_n_2;
    basis_2d_t b_m, b_n;

    real_t lambda=1.0;

    const complex_t k=2.0*pi/lambda;

    r_m_m = vector_t<real_t>(+0.3, -0.4, +0.3);
    e_m_1 = vector_t<real_t>(+0.1, +0.0, -0.7);
    e_m_2 = vector_t<real_t>(-0.3, +0.0, +0.0);
    r_m_p = vector_t<real_t>(+0.0, +0.2, +0.0);

    r_n_m = vector_t<real_t>(+0.3, -0.4, +0.3);
    e_n_1 = vector_t<real_t>(+0.1, +0.0, -0.7);
    e_n_2 = vector_t<real_t>(-0.3, +0.0, +0.0);
    r_n_p = vector_t<real_t>(+0.0, +0.2, +0.0);

    // r_n_m = vector_t<real_t>(-0.3, +0.0, +0.0);
    // e_n_1 = vector_t<real_t>(+0.1, +0.0, -0.7);
    // e_n_2 = vector_t<real_t>(+0.0, +0.2, +0.0);
    // r_n_p = vector_t<real_t>(+0.6, +0.2, +0.8);

    b_m = basis_2d_t(r_m_m, e_m_1, e_m_2, r_m_p, 1, 1);
    b_n = basis_2d_t(r_n_m, e_n_1, e_n_2, r_n_p, 1, 1);

    T.set();
    // print(phi_2d_2d(b_m, b_n, k, lambda, quadl, flag)); print(flag);
    // print(psi_2d_2d(b_m, b_n, k, lambda, quadl, flag)); print(flag);
    // print(delta_2d_2d(b_m, b_n, k, lambda, quadl, flag)); print(flag);
    print(Z_mn_2d_2d(b_m, b_n, k, eta, lambda, quadl, flag)); print(flag);
    T.unset();
    
}

void test_engine_2d_sphere_RCS(){

    // problem defintions
    const real_t GHz=1.0E+9;
    const real_t freq=0.75*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/5.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t radius=0.5;

    const size_t Ns=1001;
    complex_t E_TM, E_TE;
    real_t theta_i, phi_i;
    range_t theta_s, phi_s;
    theta_s.set(deg2rad(-180.0), deg2rad(+180.0), Ns);
    phi_s.set(deg2rad(-180.0), deg2rad(+180.0), Ns);
    theta_s.linspace();
    phi_s.linspace();

    engine_t engine;
    create_sphere(radius);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, 0, 0);

    engine.compute_Z_mn();
    // engine.load_Z_mn("data/Z_mn.bin");
    file_t file;
    sigma_t sigma;

    // theta

    E_TM = +1.0;
    E_TE = +0.0;
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_solutions();
    engine.export_currents("data/currents.pos");

    file.open("data/RCS_1.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        sigma = engine.compute_RCS(theta_s(i), phi_i);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(phi_s(i)), sigma.theta, sigma.phi);
    }
    file.close();

    // phi

    E_TM = +0.0;
    E_TE = +1.0;
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_solutions();
    
    file.open("data/RCS_2.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        sigma = engine.compute_RCS(theta_s(i), phi_i);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), sigma.theta, sigma.phi);
    }
    file.close();

    //
    theta_s.unset();
    phi_s.unset();

    engine.unset();
    
}


void test_engine_2d_sphere_near_field(){

    // problem defintions
    const real_t GHz=1.0E+9;
    const real_t freq=0.6*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/11.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t radius=0.5;

    const size_t Ns=401;
    complex_t E_TM, E_TE;
    real_t theta_i, phi_i;

    engine_t engine;
    create_sphere(radius);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, 0, 0);

    // engine.compute_Z_mn();
    engine.load_Z_mn("data/Z_mn.bin");
    file_t file;
    range_t z;
    real_t z_min, z_max, x, y;

    x = +1.0;
    y = +1.0;
    z_min = -1.0;
    z_max = +1.0;
    z.set(z_min, z_max, Ns);
    z.linspace();

    E_TM = +1.0;
    E_TE = +0.0;
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    vector_t<real_t> r;

    incident_field_t incident_field;
    vector_t<complex_t> E, H;
    file.open("data/near_field.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing near fields...");
        r = vector_t<real_t>(x, y, z(i));
        E = engine.compute_near_field_E(r);
        H = engine.compute_near_field_H(r);
        incident_field = compute_incident_field(E_TM, E_TE, theta_i, phi_i, 
            2.0*pi/lambda, sqrt(mu_0/eps_0), r);
        // E = E+incident_field.E;
        // H = H+incident_field.H;
        file.write("%21.14E ", z(i));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E ", 
            real(E.x), imag(E.x),
            real(E.y), imag(E.y),
            real(E.z), imag(E.z));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E ", 
            real(H.x), imag(H.x),
            real(H.y), imag(H.y),
            real(H.z), imag(H.z));
        file.write("\n");
    }
    file.close();

    //
    z.unset();

    engine.unset();
    
}

struct integrand_args_2d{
    triangle_t tri;
    vector_t<real_t> r, unit_vector;
};

#define EPS 1.0E-20

complex_t integrand_1_2d(const complex_t alpha, const complex_t beta, void *args_){
    assert(args_!=null);
    integrand_args_2d* args=(integrand_args_2d*)args_;
    vector_t<real_t> rho=args->tri.v[0]
                         +real(alpha)*(args->tri.v[1]-args->tri.v[0])
                         +real(beta) *(args->tri.v[2]-args->tri.v[0]);
    real_t R=mag(args->r-rho)+EPS;
    return 2.0*args->tri.area/R;
}

complex_t integrand_2_2d(const complex_t alpha, const complex_t beta, void *args_){
    assert(args_!=null);
    integrand_args_2d* args=(integrand_args_2d*)args_;
    vector_t<real_t> rho=args->tri.v[0]
                         +real(alpha)*(args->tri.v[1]-args->tri.v[0])
                         +real(beta) *(args->tri.v[2]-args->tri.v[0]);
    real_t R=mag(args->r-rho)+EPS;
    return 2.0*args->tri.area*(rho-args->r)*args->unit_vector/R;
}

complex_t integrand_3_2d(const complex_t alpha, const complex_t beta, void *args_){
    assert(args_!=null);
    integrand_args_2d* args=(integrand_args_2d*)args_;
    vector_t<real_t> rho=args->tri.v[0]
                         +real(alpha)*(args->tri.v[1]-args->tri.v[0])
                         +real(beta) *(args->tri.v[2]-args->tri.v[0]);
    real_t R=mag(args->r-rho)+EPS;
    return 2.0*args->tri.area*(1.0/pow(R, 2.0))*(args->unit_vector*unit(args->r-rho));
}

void test_engine_2d_debug(){

    //
    vector_t<real_t> v1, v2, v3;
    v1 = vector_t<real_t>(-0.5, -0.6, -0.2);
    v2 = vector_t<real_t>(+0.4, -0.6, +0.0);
    v3 = vector_t<real_t>(+0.0, +0.3, +0.2);

    triangle_t tri=triangle_t(v1, v2, v3, 0);
    vector_t<real_t> p=vector_t<real_t>(+0.2, +0.1, +1.0);

    projection_2d_para para;
    para = prjection_2d(v1, v2, v3, p);

    //
    real_t I1=0.0;
    for (size_t i=0; i<3; i++){
        real_t A=para.para_1d[i].P_0_unit*para.u[i];
        real_t B=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t D=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I1 = I1+A*(B-abs(para.d)*(C-D));
    }
    print(I1);

    vector_t<real_t> I2 = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=para.R_p[i]*para.para_1d[i].l_p-para.R_m[i]*para.para_1d[i].l_m;
        real_t B=pow(para.R_0[i], 2.0)*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        I2 = I2+0.5*(A+B)*para.u[i];
        real_t C=para.para_1d[i].P_0_unit*para.u[i];
        real_t D=para.para_1d[i].P_0*log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t E=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t F=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I2 = I2-para.d*C*(D-abs(para.d)*(E-F))*para.n;
    }
    print(I2);

    vector_t<real_t> I3 = vector_t<real_t>(0.0, 0.0, 0.0);
    for (size_t i=0; i<3; i++){
        real_t A=log((para.R_p[i]+para.para_1d[i].l_p)/(para.R_m[i]+para.para_1d[i].l_m));
        real_t B=atan2(para.para_1d[i].P_0*para.para_1d[i].l_p, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_p[i]);
        real_t C=atan2(para.para_1d[i].P_0*para.para_1d[i].l_m, 
                        pow(para.R_0[i], 2.0)+abs(para.d)*para.R_m[i]);
        I3 = I3+A*para.u[i]+sign(para.d)*(B-C)*(para.u[i]*para.para_1d[i].P_0_unit)*para.n;
    }
    print(I3);

    vector_t<real_t> x=vector_t<real_t>(1.0, 0.0, 0.0);
    vector_t<real_t> y=vector_t<real_t>(0.0, 1.0, 0.0);
    vector_t<real_t> z=vector_t<real_t>(0.0, 0.0, 1.0);

    int_t flag;
    quadl_domain_t quadl;
    quadl.set_2d(15, 1.0E-10);
    integrand_args_2d args={tri, p, p};
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), 
                                vector_t<real_t>(1.0, 0.0, 0.0), 
                                vector_t<real_t>(0.0, 1.0, 0.0)};
    
    I1 =real(quadl.integral_2d(integrand_1_2d, &args, triangle, flag)); assert(!flag);
    print(I1);

    args.unit_vector = x;
    I2.x =real(quadl.integral_2d(integrand_2_2d, &args, triangle, flag)); assert(!flag);
    args.unit_vector = y;
    I2.y =real(quadl.integral_2d(integrand_2_2d, &args, triangle, flag)); assert(!flag);
    args.unit_vector = z;
    I2.z =real(quadl.integral_2d(integrand_2_2d, &args, triangle, flag)); assert(!flag);
    print(I2);

    args.unit_vector = x;
    I3.x =real(quadl.integral_2d(integrand_3_2d, &args, triangle, flag)); assert(!flag);
    args.unit_vector = y;
    I3.y =real(quadl.integral_2d(integrand_3_2d, &args, triangle, flag)); assert(!flag);
    args.unit_vector = z;
    I3.z =real(quadl.integral_2d(integrand_3_2d, &args, triangle, flag)); assert(!flag);
    print(I3);

}

void test_engine_2d_sheet_near_field(){

    // problem defintions
    const real_t GHz=1.0E+9;
    const real_t freq=2.0*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/11.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t Lx=0.6, Ly=0.4;

    const size_t Ns=201;
    complex_t E_TM, E_TE;
    real_t theta_i, phi_i;

    engine_t engine;
    create_sheet(Lx, Ly, clmax);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, 0, 0);

    engine.compute_Z_mn();
    // engine.load_Z_mn("data/Z_mn.bin");
    file_t file;
    range_t z;
    real_t z_min, z_max, x, y;

    x = +1.0;
    y = +1.0;
    z_min = -1.0;
    z_max = +1.0;
    z.set(z_min, z_max, Ns);
    z.linspace();

    E_TM = +1.0;
    E_TE = +0.0;
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_currents("data/currents.pos");
    vector_t<real_t> r;

    vector_t<complex_t> E, H;
    file.open("data/near_field.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing near fields...");
        r = vector_t<real_t>(x, y, z(i));
        E = engine.compute_near_field_E(r);
        H = engine.compute_near_field_H(r);
        file.write("%21.14E ", z(i));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E ", 
            real(E.x), imag(E.x),
            real(E.y), imag(E.y),
            real(E.z), imag(E.z));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E ", 
            real(H.x), imag(H.x),
            real(H.y), imag(H.y),
            real(H.z), imag(H.z));
        file.write("\n");
    }
    file.close();

    //
    z.unset();

    engine.unset();
    
}

void test_engine_2d_box_RCS(){

    // problem defintions
    const real_t mm=1.0E-3;
    const real_t GHz=1.0E+9;
    const real_t freq=0.43*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/5.0;
    const complex_t mu_b=1.0, eps_b=1.0;

    const size_t Ns=201;
    complex_t E_TM, E_TE;
    real_t theta_i, phi_i;
    range_t theta_s;
    theta_s.set(deg2rad(-180.0), deg2rad(+180.0), Ns);
    theta_s.linspace();

    engine_t engine;
    create_box();
    engine.set(freq, mu_b, eps_b, clmax/mm, mm, 0, 0);

    engine.compute_Z_mn();
    // engine.load_Z_mn("data/Z_mn.bin");
    file_t file;
    sigma_t sigma;

    // theta

    E_TM = +1.0;
    E_TE = +0.0;

    file.open("data/RCS_1.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "RCS VV...");
        theta_i = theta_s(i);
        phi_i = deg2rad(0.0);
        engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
        engine.compute_I_n();
        sigma = engine.compute_RCS(theta_s(i), phi_i);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_i), sigma.theta, sigma.phi);
    }
    file.close();

    // phi

    E_TM = +0.0;
    E_TE = +1.0;
    
    file.open("data/RCS_2.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "RCS HH...");
        theta_i = theta_s(i);
        phi_i = deg2rad(0.0);
        engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
        engine.compute_I_n();
        sigma = engine.compute_RCS(theta_s(i), phi_i);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_i), sigma.theta, sigma.phi);
    }
    file.close();

    //
    theta_s.unset();

    engine.unset();
    
}

void test_engine_2d_box_near_field(){

    // problem defintions
    const real_t mm=1.0E-3;
    const real_t GHz=1.0E+9;
    const real_t freq=4.0*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/11.0;
    const complex_t mu_b=1.0, eps_b=1.0;

    const size_t Ns=401;
    complex_t E_TM, E_TE;
    real_t theta_i, phi_i;
    
    engine_t engine;
    create_box();
    engine.set(freq, mu_b, eps_b, clmax/mm, mm, 0, 0);

    engine.compute_Z_mn();
    // engine.load_Z_mn("data/Z_mn.bin");
    file_t file;
    range_t z;
    real_t z_min, z_max, x, y;

    x = +100.0*mm;
    y = +100.0*mm;
    z_min = -50.0*mm;
    z_max = +50.0*mm;
    z.set(z_min, z_max, Ns);
    z.linspace();

    E_TM = +1.0;
    E_TE = +0.0;
    theta_i = deg2rad(+0.0);
    phi_i = deg2rad(+0.0);
    engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_currents("data/currents.pos");
    vector_t<real_t> r;

    vector_t<complex_t> E, H;
    file.open("data/near_field.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing near fields...");
        r = vector_t<real_t>(x, y, z(i));
        E = engine.compute_near_field_E(r);
        H = engine.compute_near_field_H(r);
        file.write("%21.14E ", z(i));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E ", 
            real(E.x), imag(E.x),
            real(E.y), imag(E.y),
            real(E.z), imag(E.z));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E ", 
            real(H.x), imag(H.x),
            real(H.y), imag(H.y),
            real(H.z), imag(H.z));
        file.write("\n");
    }
    file.close();

    //
    z.unset();

    engine.unset();
    
}


void test_engine_2d_sphere_near_field_2d(){

    // problem defintions
    const real_t GHz=1.0E+9;
    const real_t freq=0.5*GHz;
    const real_t lambda=c_0/freq;
    const real_t clmax=lambda/11.0;
    const complex_t mu_b=1.0, eps_b=1.0;
    const real_t radius=0.5;

    const size_t Ns_x=401, Ns_z=401;
    complex_t E_TM, E_TE;
    real_t theta_i, phi_i;

    engine_t engine;
    create_sphere(radius);
    engine.set(freq, mu_b, eps_b, clmax, 1.0, 0, 0);

    engine.compute_Z_mn();
    // engine.load_Z_mn("data/Z_mn.bin");
    file_t file_x, file_y, file_data;
    range_t x, z;
    real_t z_min, z_max, x_min, x_max, y;
    const real_t range=1.0;

    x_min = -range;
    x_max = +range;
    y = +0.0;
    z_min = -range;
    z_max = +range;
    x.set(x_min, x_max, Ns_x);
    z.set(z_min, z_max, Ns_z);
    x.linspace();
    z.linspace();

    E_TM = +1.0;
    E_TE = +0.0;
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    engine.compute_V_m_incident(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    vector_t<real_t> r;

    vector_t<complex_t> E, H;
    incident_field_t incident_field;
    file_x.open("data/near_field_2d_x.txt", 'w');
    file_y.open("data/near_field_2d_y.txt", 'w');
    file_data.open("data/near_field_2d_data.txt", 'w');

    const size_t max_line=200;
    char *msg=(char*)calloc(max_line, sizeof(char));
    size_t counter=0;
    for (size_t i=0; i<Ns_x; i++){
        for (size_t j=0; j<Ns_z; j++){
            sprintf(msg, "computing near fields (%zu, %zu)...", i, j);
            progress_bar(counter, Ns_x*Ns_z, msg);
            counter++;
            r = vector_t<real_t>(x(i), y, z(j));
            incident_field = compute_incident_field(E_TM, E_TE, theta_i, phi_i, 
                                                    2.0*pi/lambda, sqrt(mu_0/eps_0), r);
            E = engine.compute_near_field_E(r);
            H = engine.compute_near_field_H(r);
            E = E+incident_field.E;
            H = H+incident_field.H;
            file_x.write("%21.14E ", x(i));
            file_y.write("%21.14E ", z(j));
            vector_t<complex_t> E_field=vector_t<complex_t>(E.x, E.y, E.z);
            vector_t<complex_t> H_field=vector_t<complex_t>(H.x, H.y, H.z);
            // file_data.write("%21.14E ", sqrt(pow(real(E.x), 2.0)
            //                                 +pow(real(E.y), 2.0)
            //                                 +pow(real(E.z), 2.0)));
            file_data.write("%21.14E ", sqrt(pow(real(H.x), 2.0)
                                            +pow(real(H.y), 2.0)
                                            +pow(real(H.z), 2.0)));
        }
        file_x.write("\n");
        file_y.write("\n");
        file_data.write("\n");
    }
    free(msg);
    file_x.close();
    file_y.close();
    file_data.close();

    //
    x.unset();
    z.unset();

    engine.unset();
    
}