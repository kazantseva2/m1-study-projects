#ifndef CLASSES_H
#define CLASSES_H

#include <vector>
#include <string>
#include <cmath> 
#include <random> 
#include <algorithm>
#include <iomanip> 


using namespace std;

// Перечисление металлов
enum Metal {
    ARGENTUM=0,
    FERRUM=1
};

// Класс для работы с векторами
struct Vector {
    double x, y, z;

    Vector(double x0, double y0, double z0);
    Vector();

    void print() const;
    void modf();
    bool equals(const Vector& v) const;
    double length() const;
    Vector add(double a, double b, double c) const;
    Vector multiplyByVector(const Vector& coefs) const;
    Vector multiplyByMatrix(const vector<vector<double>>& matrix) const;
    bool operator==(const Vector& other) const;
    Vector operator-(const Vector& other) const;
};

// Класс, представляющий атом
struct Atom {
    Vector coord;
    Metal type;

    Atom(Vector v, Metal metal);
};

// Класс, представляющий решётку
struct Lattice {
    int N, K; // количество ГЦК в решетке и общее количество атомов в решетке
    Vector a;
    double cutoff, V;
    vector<Atom> atoms;
    vector<Vector> neighbours;
    vector<Vector> fccBasis = {
        Vector(0.0, 0.0, 0.0),
        Vector(0.5, 0.5, 0.0),
        Vector(0.5, 0.0, 0.5),
        Vector(0.0, 0.5, 0.5)
    };

    Lattice(int N_init, double a_init, double cutoff_init = 1.7);
    Lattice(const Lattice& original, vector<vector<double>> deltaMatrix);
    
    bool isAtomInFCC(Vector atom);
    
    vector<Vector> initNeighbours();
    vector<Atom> initLattice(int N, Metal metal);
    int getIndexAtom(Vector v) const;
};

enum Structure {
    BASE=0,
    EXTRA=1,
    SOLUTION=2
};

struct Parameters {
    double xi = 1.178;
    double r0 = 4.085 / sqrt(2);
    double q = 3.139;
    double A1 = 0, A0 = 0.1028, pr = 10.928;

    double error;

    void print() const;
    void operator*=(double scalar);
    void operator+=(Parameters other);
    Parameters operator*(double scalar) const;
    Parameters operator-(Parameters other) const;
    Parameters operator+(Parameters other) const;
    void randInit();
};

struct ParametersGroup {
    Parameters paramAA;
    Parameters paramAB;
    Parameters paramBB;
};

// Класс для расчёта энергии
class Energy {
public:
    ParametersGroup pg;
    double computeBindingEnergy(double r, Structure s) const;
    double computeRepulsionEnergy(double r, Structure s) const;


    Energy();
    Energy(ParametersGroup &paramGroup);
    double computeEnergy(const Lattice& lattice, bool vacuum = false) const;
    double computeCohesionEnergy(const Lattice& lt) const;
    double computeSolubilityEnergy(Lattice& lt);
    double computeEnergyOfDimerInSurface(Lattice& lt);
    double computeEnergyOfDimerOnSurface(Lattice& lt);
    double computeEnergyOne(double r, Structure s) const;
};

struct ElasticityValues {
    double B;  
    double C11;     
    double C12;     
    double C44;     
};

class Elasticity {
    double derivative(Energy &energy, Lattice &originalLattice, double delta, vector<vector<double>> &deformationMatrixMinus, vector<vector<double>> &deformationMatrixPlus) const;
    double BDerivative(Energy &energy, Lattice &originalLattice, double delta) const;
    double C11Derivative(Energy &energy, Lattice &originalLattice, double delta) const;
    double C12Derivative(Energy &energy, Lattice &originalLattice, double delta) const;
    double C44Derivative(Energy &energy, Lattice &originalLattice, double delta) const;

public:
    ElasticityValues getElasticity(Energy &energy, Lattice &lt, double delta) const;
    void elasticityPrint(Energy &energy, Lattice &lt, double delta) const;
};



class Optimize {
public:
    double a_opt = 4.085;
    double E_coh_opt = -2.96;
    double B_opt = 1.08;
    double C11_opt = 1.32;
    double C12_opt = 0.97;
    double C44_opt = 0.51;
    double E_sol_opt = 0.974;
    double E_in_dim_opt = -0.01;
    double E_on_dim_opt = -0.72;
    
    ParametersGroup pg_opt;

    double alpha = 1.0, beta = 0.5, gamma = 2, sigma = 0.5;

    vector<Parameters> simplex;

    double mse(double f_val, double f_opt, double w);
    void initializeSimplex(Lattice &lt, Structure st, int n = 6);
    void sortSimplexByValue();
    Parameters getCenter();
    double errorAA(Parameters &p, Lattice &lt);
    double errorAB(Parameters &p, Lattice &lt);
    double errorBB(Parameters &p, Lattice &lt);
    void errorFun(Parameters &p, Lattice &lt, Structure st);

    
    ParametersGroup optimize(Lattice &lt, int n, Structure st,  double eps);
    double printResult(Lattice &lt);
};

#endif // CLASSES_H
