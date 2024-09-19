#ifndef __ENGINE_HPP__
#define __ENGINE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
#include "shape.hpp"
#include "projection.hpp"
//
#include "shared_definitions.hpp"
#include "R_mn.hpp"
#include "integrand_projection_1d.hpp"
// #include "integrand_projection_2d.hpp"
// #include "integrand_projection_3d.hpp"
#include "psi_integrand.hpp"
#include "phi_integrand.hpp"
#include "delta_integrand.hpp"
#include "Z_mn.hpp"
#include "incident_field.hpp"
#include "scattered_field.hpp"

// Definitions
#define max_system_size 20000

struct port_t{
    size_t index=0;
    complex_t V=0.0;
    complex_t Z=0.0;
    int_t pg=0;
    real_t L=0.0, W=0.0;
    vector_t<real_t> p=vector_t<real_t>(0.0, 0.0, 0.0);
    port_t(){}
};

struct sigma_t{
    real_t theta=0.0, phi=0.0;
};

struct far_field_t{
    complex_t theta=0.0, phi=0.0;
};

class engine_t{
    private:
        const size_t max_line_length=200;
        //
        int_t is_engine_set=false;
        int_t is_Z_mn_allocated=false;
        int_t is_V_m_allocated=false;
        int_t is_I_n_allocated=false;
        int_t is_Z_mn_calculated=false;
        int_t is_V_m_calculated=false;
        int_t is_I_n_calculated=false;
        int_t is_port_list_allocated=false;
        //
        shape_t shape;
        matrix_t<complex_t> Z_mn, V_m, I_n;
        real_t freq=0.0, lambda=0.0;
        complex_t mu_b=1.0, eps_b=1.0;
        complex_t k_b=0.0, eta_b=0.0;
        real_t a=0.0;
        size_t N_basis_1d=0;
        size_t N_basis_2d=0;
        size_t N_basis_3d=0;
        size_t N=0;
        //
        quadl_domain_t quadl;
        const size_t k_max_1d=15; const real_t tol_1d=1.0E-4;
        const size_t k_max_2d=15; const real_t tol_2d=1.0E-4;
        const size_t k_max_3d=15; const real_t tol_3d=1.0E-4;
        //
        size_t N_ports=0;
        port_t *port_list=null;
    public:        
        void save_Z_mn(const char *filename);
        void load_Z_mn(const char *filename);
        void get_basis_info(size_t &N_1d, size_t &N_2d, size_t &N_3d){
            N_1d = this->N_basis_1d;
            N_2d = this->N_basis_2d;
            N_3d = this->N_basis_3d;
        }
        engine_t(){}
        ~engine_t(){}
        void set(const real_t freq, const complex_t mu_b, const complex_t eps_b, 
            const real_t clmax, const real_t unit_metric, const real_t a, const size_t N_ports);
        void compute_Z_mn();
        void unset();
        void export_solutions();
        void assign_port(const size_t index, const complex_t V, const complex_t Z, const int_t pg, 
            const vector_t<real_t> p, const real_t L, const real_t W);
        void compute_V_m_ports();
        void compute_V_m_incident(const complex_t E_TM, const complex_t E_TE, const real_t theta_i, const real_t phi_i);
        void compute_I_n();
        complex_t compute_Z_in(const size_t port_index);
        complex_t compute_S_mutual(const size_t port_index);
        void compute_S_matrix(matrix_t<complex_t> &S_matrix, const complex_t Z_0);
        //
        sigma_t compute_RCS(const real_t theta_i, const real_t phi_i);
        far_field_t compute_far_field(const real_t theta_i, const real_t phi_i);
        vector_t<complex_t> compute_near_field_E(const vector_t<real_t> r);
        vector_t<complex_t> compute_near_field_H(const vector_t<real_t> r);
        void export_currents(const char *filename);
        void set_material(const int_t pg, const complex_t mu, const complex_t eps){
            this->shape.set_material(pg, mu, eps);
        }
};

// Functions

#endif
