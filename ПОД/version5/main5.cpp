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

        double secondDerivative = (E1 - 2*E2 + E3)/(delta*delta);
        return secondDerivative;
    }

    double EC11(Energy &energy, Lattice &originalLattice, double delta) const {
        return energy.computeEnergy(Lattice(originalLattice, Vector(1+delta, 1+delta, 1)));
    }

    double EC12(Energy &energy, Lattice &originalLattice, double delta) const {
        return energy.computeEnergy(Lattice(originalLattice, Vector(1+delta, 1-delta, 1)));
    }

    void elasticityPrint(Energy &energy, Lattice &lt, double delta) const {
        double BDeriv = BDerivative(energy, lt, delta);
        double EC11value = EC11(energy, lt, delta);
        double EC12value = EC12(energy, lt, delta);

        double Bvalue = (2*BDeriv)/(9*lt.V);
        double C11 = (EC11value + EC12value)/(2 * lt.V *delta*delta);
        double C12 = (EC11value - EC12value)/(2 * lt.V*delta*delta);

        cout << "BDeriv = " << BDeriv << endl;
        cout << "EC11value = " << EC11value << endl;
        cout << "EC12value = " << EC12value << endl;

        cout << "B = " << Bvalue << endl;
        cout << "C11 = " << C11 << endl;
        cout << "C12 = " << C12 << endl;
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

    return 0;
}