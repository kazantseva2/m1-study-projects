#include <iostream>
#include <cmath>  
#include <fstream>
#include <vector>

using namespace std;

enum Metal {
    ARGENTUM,
    FERRUM
};

struct Vector {
    double x, y, z;

    Vector(double x0, double y0, double z0) {
        x = x0;
        y = y0;
        z = z0;
    }

    Vector() {}


    void print() const {
        cout << "(" << x << ", " << y  << ", "<< z <<")";
    }

    void modf() {
        x = abs(x - static_cast<int>(x));
        y = abs(y - static_cast<int>(y));
        z = abs(z - static_cast<int>(z));
    }

    bool equals(const Vector& v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    double length() const {
        return sqrt(x*x + y*y + z*z);
    }

    Vector add(double a, double b, double c) const {
        return Vector(x + a, y + b, z + c);
    }

    // Умножение на вектор
    Vector multiplyByVector(const Vector& coefs) const {
        return Vector(x * coefs.x, y * coefs.y, z * coefs.z);
    }

    // Умножение на матрицу 3*3
    Vector multiplyByMatrix(const vector<vector<double>>& matrix) const {
        return Vector(
            x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2],
            x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2],
            x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2]
        );
    }
};

struct Atom {
    Vector coord;
    Metal type;

    Atom(Vector v, Metal metal) {
        coord = v;
        type = metal;
    }
};

struct Lattice {
    int N;
    int K;
    Vector a;
    double cutoff;
    double V;
    vector<Atom> atoms;
    vector<Vector> neighbours;
    vector<Vector> fccBasis = {
        Vector(0.0, 0.0, 0.0),
        Vector(0.5, 0.5, 0.0),
        Vector(0.5, 0.0, 0.5),
        Vector(0.0, 0.5, 0.5)
    };

    // Стандартный конструктор
    Lattice(int N_init, double a_init, double cutoff_init = 1.7) {
        N = N_init;
        K = 4*N*N*N;
        a = Vector(a_init, a_init, a_init);
        V = a_init * a_init * a_init * N * N * N;
        cutoff = cutoff_init;
        neighbours = initNeighbours(cutoff);
        atoms = initLattice(N, Metal::ARGENTUM);
    }

    // Конструктор для деформированной матрицы (есть копирование атомов и соседей) + (не регулируем cutoff)
    Lattice(const Lattice& original, Vector delta) {
        N = original.N;
        K = original.K;
        a = original.a.multiplyByVector(delta);
        V = a.x * a.y * a.z * N * N * N;
        cutoff = original.cutoff;
        fccBasis = original.fccBasis;
        neighbours = initNeighbours(cutoff);
        atoms = initLattice(N, Metal::ARGENTUM);
    }


    bool isAtomInFCC(Vector atom) {
        atom.modf();
        for (const Vector& basisAtom : fccBasis) {
            if (atom.equals(basisAtom)) {
                return true;
            }
        }

        return false;
    }

    vector<Vector> initNeighbours(double cutoff, double delta = 0) {
        vector<Vector> neighbours;
        double coordCutoff = static_cast<int>(cutoff / 0.5) * 0.5;
        Vector center = Vector(0, 0, 0);
        Vector nb = Vector(center);
        int n_around = 0;
        for (nb.x = center.x - coordCutoff; nb.x <= center.x + cutoff; nb.x += 0.5) { 
            for (nb.y = center.y - coordCutoff; nb.y <= center.y + cutoff; nb.y += 0.5) {
                for (nb.z = center.z - coordCutoff; nb.z <= center.z + cutoff; nb.z += 0.5) {
                    if (!isAtomInFCC(nb) || center.equals(nb) || nb.length() > cutoff) {
                        continue;
                    }
                    neighbours.push_back(Vector(nb));
                }
            }
        }
        return neighbours;
    }

    vector<Atom> initLattice(int N, Metal metal) {
        vector<Atom> atoms;
        for (Vector basisAtom : fccBasis) {
            for (int k = 0; k < N; k++) { 
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < N; i++) {
                        Vector atom_coord = basisAtom.add(i, j, k);
                        atoms.push_back(Atom(atom_coord, metal));
                    }
                }
            }
        }
        return atoms;   
    }
};

class Energy {
    double xi = 1.178;
    double r0 = 4.085 / sqrt(2);
    double q = 3.139;
    double A1 = 0, A0 = 0.1028, p = 10.928;

    double computeBindingEnergy(double r) const {
        double tmp = r/r0 - 1;
        return exp(-2 * q * tmp);
    }

    double computeRepulsionEnergy(double r) const {
        double tmp = r/r0 - 1;
        return (A1 * (r - r0) + A0) * exp(-1 * p * tmp);
    }

    public:
    double computeEnergy(const Lattice& lt) const {
        double Eb, Er, E = 0;
        for (const Atom atom : lt.atoms) {
            Eb = Er = 0;
            for (const Vector nb : lt.neighbours) { 
                double r = nb.multiplyByVector(lt.a).length();
                if (nb.x == 0 && nb.y == 0) {
                   // nb.print();
                    // cout << " r = " << r <<endl;
                }
                Eb += computeBindingEnergy(r);
                Er += computeRepulsionEnergy(r);
            }
            E += Er - xi*sqrt(Eb);
        }
        return E;
    }

    double computeCohesionEnergy(const Lattice& lt) const {
        return computeEnergy(lt)/lt.K;
    }
};

struct Elasticity {
    double BDerivative(Energy &energy, Lattice &originalLattice, double delta) const {
        // vector<vector<double>> deformationMatrixMinus = {
        //     {1 - delta, 0.0, 0.0},
        //     {0.0, 1 - delta, 0.0},
        //     {0.0, 0.0, 1 - delta}
        // };

        // vector<vector<double>> deformationMatrixPlus = {
        //     {1 + delta, 0.0, 0.0},
        //     {0.0, 1 + delta, 0.0},
        //     {0.0, 0.0, 1 + delta}
        // };

        double E1 = energy.computeEnergy(Lattice(originalLattice, Vector(1-delta, 1-delta, 1-delta)));
        double E2 = energy.computeEnergy(originalLattice);
        double E3 = energy.computeEnergy(Lattice(originalLattice, Vector(1+delta, 1+delta, 1+delta)));
        cout << "BDerivative: E- = " <<  E1  << ", ";
        cout << "BDerivative: E0 = " <<  E2  << ", ";
        cout << "BDerivative: E+ = " <<  E3  << endl;

        double secondDerivative = (E1 - 2*E2 + E3)/(delta*delta);
        return secondDerivative;
    }

    double C11Derivative(Energy &energy, Lattice &originalLattice, double delta) const {
        double E1 = energy.computeEnergy(Lattice(originalLattice, Vector(1-delta, 1-delta, 1)));
        double E2 = energy.computeEnergy(originalLattice);
        double E3 = energy.computeEnergy(Lattice(originalLattice, Vector(1+delta, 1+delta, 1)));
        cout << "C11E = " << E3<<endl;

        double secondDerivative = (E1 - 2*E2 + E3)/(delta*delta);
        return secondDerivative;
    }

    double C12Derivative(Energy &energy, Lattice &originalLattice, double delta) const {
        double E1 = energy.computeEnergy(Lattice(originalLattice, Vector(1-delta, 1+delta, 1)));
        double E2 = energy.computeEnergy(originalLattice);
        double E3 = energy.computeEnergy(Lattice(originalLattice, Vector(1+delta, 1-delta, 1)));
        cout << "C12E = " << E3<<endl;

        double secondDerivative = (E1 - 2*E2 + E3)/(delta*delta);
        return secondDerivative;
    }

    void elasticityPrint(Energy &energy, Lattice &lt, double delta) const {
        double C11Deriv = C11Derivative(energy, lt, delta);
        double C12Deriv = C12Derivative(energy, lt, delta);
        double BDeriv = BDerivative(energy, lt, delta);

        double C11 = (C11Deriv + C12Deriv)/(lt.V * 2);
        double C12 = (C11Deriv - C12Deriv)/(lt.V * 2);
        double Bvalue = (2*BDeriv)/(9*lt.V);

        cout << "BDeriv = " << BDeriv << endl;
        cout << "C11Deriv = " << C11Deriv << endl;
        cout << "C12Deriv = " << C12Deriv << endl;

        cout << "B = " << Bvalue << endl;
        cout << "C11 = " << C11 << endl;
        cout << "C12 = " << C12 << endl;

        double d = 0.624*2;
        cout << "B = " << Bvalue/d << endl;
        cout << "C11 = " << C11/d << endl;
        cout << "C12 = " << C12/d << endl;
    }

};

int main()
{
    int N = 1;
    Lattice lattice = Lattice(N, 4.085);
    Energy energy = Energy();
    cout << "Энергия когезии: " << energy.computeCohesionEnergy(lattice) << endl;

    int i =0;
    for (auto nb : lattice.neighbours) {
        i++;
    }
    cout << "Количество соседей, которое учитывается: " << i << endl;

    Elasticity elasticity = Elasticity();
    cout << "Модуль упругости: "<< endl;
    elasticity.elasticityPrint(energy, lattice, 0.01);


    Lattice lattice1 = Lattice(N, 4.04415);
    Lattice lattice0 = Lattice(N, 4.085);
    Lattice lattice2 = Lattice(N, 4.12585);
    double V = 4.085*4.085*4.085;
    double delta = 0.01;
    double E1, E0, E2;
    E1 = energy.computeEnergy(lattice1);
    E0 = energy.computeEnergy(lattice0);
    E2 = energy.computeEnergy(lattice2);
    cout << "E- : " <<  E1  << endl;
    cout << "E0 : " <<  E0  << endl;
    cout << "E+ : " <<  E2  << endl;
    cout << "cohE- : " <<  energy.computeCohesionEnergy(lattice1)  << endl;
    cout << "cohE0 : " <<  energy.computeCohesionEnergy(lattice0)  << endl;
    cout << "cohE+ : " <<  energy.computeCohesionEnergy(lattice2)  << endl;
    double res = 2*(E1 - 2*E0 + E2)/(9*V *delta*delta);
    cout << "res = " << res << endl;

    return 0;
}