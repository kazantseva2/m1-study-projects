#include "classes.h"
#include <iostream>


using namespace std;

double Elasticity::derivative(Energy &energy, Lattice &originalLattice, double delta, vector<vector<double>> &deformationMatrixMinus, vector<vector<double>> &deformationMatrixPlus) const {
    double E1 = energy.computeEnergy(Lattice(originalLattice, deformationMatrixMinus));
    double E2 = energy.computeEnergy(originalLattice);
    double E3 = energy.computeEnergy(Lattice(originalLattice, deformationMatrixPlus));

    double secondDerivative = (E1 - 2*E2 + E3)/(delta*delta);
    return secondDerivative;
}

double Elasticity::BDerivative(Energy &energy, Lattice &originalLattice, double delta) const {
    vector<vector<double>> deformationMatrixMinus = {
        {1 - delta, 0.0, 0.0},
        {0.0, 1 - delta, 0.0},
        {0.0, 0.0, 1 - delta}
    };

    vector<vector<double>> deformationMatrixPlus = {
        {1 + delta, 0.0, 0.0},
        {0.0, 1 + delta, 0.0},
        {0.0, 0.0, 1 + delta}
    };

    return derivative(energy, originalLattice, delta, deformationMatrixMinus, deformationMatrixPlus);
}

double Elasticity::C11Derivative(Energy &energy, Lattice &originalLattice, double delta) const {
    vector<vector<double>> deformationMatrixMinus = {
        {1 - delta, 0.0, 0.0},
        {0.0, 1 - delta, 0.0},
        {0.0, 0.0, 1.0}
    };

    vector<vector<double>> deformationMatrixPlus = {
        {1 + delta, 0.0, 0.0},
        {0.0, 1 + delta, 0.0},
        {0.0, 0.0, 1.0}
    };

    return derivative(energy, originalLattice, delta, deformationMatrixMinus, deformationMatrixPlus);
}

double Elasticity::C12Derivative(Energy &energy, Lattice &originalLattice, double delta) const {
    vector<vector<double>> deformationMatrixMinus = {
        {1 - delta, 0.0, 0.0},
        {0.0, 1 + delta, 0.0},
        {0.0, 0.0, 1.0}
    };

    vector<vector<double>> deformationMatrixPlus = {
        {1 + delta, 0.0, 0.0},
        {0.0, 1 - delta, 0.0},
        {0.0, 0.0, 1.0}
    };

    return derivative(energy, originalLattice, delta, deformationMatrixMinus, deformationMatrixPlus);
}

double Elasticity::C44Derivative(Energy &energy, Lattice &originalLattice, double delta) const {
    vector<vector<double>> deformationMatrixMinus = {
        {1.0, 0-delta, 0.0},
        {0-delta, 1.0, 0.0},
        {0.0, 0.0, 1/(1-delta*delta)}
    };

    vector<vector<double>> deformationMatrixPlus = {
        {1.0, delta, 0.0},
        {delta, 1.0, 0.0},
        {0.0, 0.0, 1/(1-delta*delta)}
    };

    return derivative(energy, originalLattice, delta, deformationMatrixMinus, deformationMatrixPlus);
}

ElasticityValues Elasticity::getElasticity(Energy &energy, Lattice &lt, double delta) const {
    double C11Deriv = C11Derivative(energy, lt, delta);
    double C12Deriv = C12Derivative(energy, lt, delta);
    double BDeriv = BDerivative(energy, lt, delta);
    double C44Deriv = C44Derivative(energy, lt, delta);

    double C11 = (C11Deriv + C12Deriv)/(lt.V * 2);
    double C12 = (C11Deriv - C12Deriv)/(lt.V * 2);
    double Bvalue = (2*BDeriv)/(9*lt.V);
    double C44 = C44Deriv/(lt.V * 8);

    double d = 0.624*2;

    return {Bvalue/d, C11/d, C12/d, C44/d};
}

void Elasticity::elasticityPrint(Energy &energy, Lattice &lt, double delta) const {
    ElasticityValues ElasticityValues = getElasticity(energy, lt, delta);

    double C11 = ElasticityValues.C11;
    double C12 = ElasticityValues.C12;
    double Bvalue = ElasticityValues.B;
    double C44 = ElasticityValues.C44;


    cout << "B = " << Bvalue << endl;
    cout << "C11 = " << C11 << endl;
    cout << "C12 = " << C12 << endl;
    cout << "C44 = " << C44 << endl;
}