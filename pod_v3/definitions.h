#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <vector>
#include <cmath> 
#include <algorithm>
#include <random> 
#include <omp.h>

using namespace std;


struct Vector {
    double x, y, z;

    Vector(double x0, double y0, double z0);
    Vector();

    Vector add(double a, double b, double c) const;
    bool operator==(Vector other) const;
    Vector operator-(Vector other) const;
    Vector operator*(Vector other) const;
    double length() const;
    double mod_x(double x, double n) const;
    Vector mod(double n, bool vacuum) const;
    Vector operator*(const vector<vector<double>>& matrix) const;
};

enum Metal {
    ARGENTUM=0,
    VANADIUM=1
};

struct Atom {
    Vector coord;
    Metal type;

    Atom(Vector v, Metal metal);
    bool operator==(Atom other) const;
};

struct Lattice {
    int N, K;
    Vector a;
    double cutoff, V;
    vector<Vector> fcc_basis = {
        Vector(0.0, 0.0, 0.0),
        Vector(0.5, 0.5, 0.0),
        Vector(0.5, 0.0, 0.5),
        Vector(0.0, 0.5, 0.5)
    };
    vector<Atom> atoms;

    Lattice(int N_init, double a_init, double cutoff_init = 1.7);
    
    vector<Atom> init_lattice(int N, Metal metal);
    int get_index(const Vector& coord) const;
};

enum Interaction {
    BASE=0,
    EXTRA=1,
    SOLUTION=2
};

struct Parameters {
    double xi = 1.178;
    double r0 = 4.085 / sqrt(2);
    double q = 3.139;
    double A1 = 0, A0 = 0.1028, pr = 10.928;

    double error_value;

    void print() const;
    void rand_init();
    void operator*=(double scalar);
    void operator+=(Parameters other);
    Parameters operator*(double scalar) const;
    Parameters operator-(Parameters other) const;
    Parameters operator+(Parameters other) const;
};

struct ParametersGroup {
    Parameters paramAA;
    Parameters paramAB;
    Parameters paramBB;
};

struct ElasticityValues {
    double B;  
    double C_11;     
    double C_12;     
    double C_44;     
};

struct Calculator {
    ParametersGroup& params;
    Lattice& lattice;
    double delta;

    Calculator(Lattice& lattice_init, ParametersGroup& pg_init, double delta_init);

    Interaction get_interaction(const Atom& atom1, const Atom& atom2) const;
    Vector get_mod_dist_vector(const Atom& atom1, const Atom& atom2, bool vacuum) const;
    double compute_binding_energy(double r, Interaction interaction) const;
    double compute_repulsion_energy(double r, Interaction interaction) const;
    double compute_energy(bool vacuum = false) const;
    double compute_energy(vector<vector<double>> &deformation_matrix, bool vacuum = false) const;
    double compute_cohesion_energy() const;

    double derivative(vector<vector<double>> &deformation_matrix_minus, vector<vector<double>> &deformation_matrix_plus);
    double derivative_B();
    double derivative_C_11();
    double derivative_C_12();
    double derivative_C_44();
    ElasticityValues get_elasticity();
    void elasticity_print();

    double compute_solubility_energy();
    double compute_energy_of_dimer_in_surface();
    double compute_energy_of_dimer_on_surface();

};

struct Optimize {
    double a_opt = 4.085;
    double E_coh_opt = -2.96;
    double B_opt = 1.08;
    double C_11_opt = 1.32;
    double C_12_opt = 0.97;
    double C_44_opt = 0.51;
    double E_sol_opt = 0.497;
    double E_in_dim_opt = 0.22;
    double E_on_dim_opt = -0.36;
    
    Calculator& calculator;

    double alpha = 1.0, beta = 0.5, gamma = 2, sigma = 0.5;

    vector<Parameters> simplex;

    Optimize(Calculator& calculator_init);

    double mse(double f_val, double f_opt, double w);
    double error_AA(Parameters &p);
    double error_AB(Parameters &p);
    double error_BB(Parameters &p);
    void error_fun(Parameters &p, Interaction interact);
    void sort_simplex_by_value();
    Parameters get_center();
    void initialize_simplex(Interaction interact, int n);
    void optimize(Interaction interact, int n, double eps);
    void print_parameters() const;
    void print_result() const;
};


#endif // CLASSES_H