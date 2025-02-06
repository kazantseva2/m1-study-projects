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

        void print() {
        int N = 3;
        double a = x, b = y, c = z;
        if (a > N/2) {
            a -=N;
        }
        if (b > N/2) {
            b -=N;
        }
        if (c > N/2) {
            c -=N;
        }
        cout << "(" << a << ", " << b  << ", "<< c <<")";
    }

    // void print() {
    //     cout << "(" << x << ", " << y  << ", "<< z <<")";
    // }

    void modf() {
        x = abs(x - static_cast<int>(x));
        y = abs(y - static_cast<int>(y));
        z = abs(z - static_cast<int>(z));
    }

    bool equals(const Vector& v) {
        return x == v.x && y == v.y && z == v.z;
    }

    // double length() {
    //     return sqrt(x*x + y*y + z*z);
    // }

    Vector add(double a, double b, double c){
        return Vector(x + a, y + b, z + c);
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
    double a;
    double cutoff = 1.7;
    vector<Atom> atoms; // хранятся в относительных координат
    vector<Vector> fccBasis = {
        Vector(0.0, 0.0, 0.0),
        Vector(0.5, 0.5, 0.0),
        Vector(0.5, 0.0, 0.5),
        Vector(0.0, 0.5, 0.5)
    };

    Lattice(int N_init, double a_init) {
        N = N_init;
        a = a_init;
        atoms = initLattice(N, Metal::ARGENTUM);
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

    double cycle_abs(double x) {
        x = abs(x);
        if (x > N/2) {
            return N - x;
        }
        return x;
    }

    double dist(Atom& atom1, Atom& atom2) {
        Vector coord1 = atom1.coord;
        Vector coord2 = atom2.coord;
        double x = cycle_abs(coord1.x - coord2.x);
        double y = cycle_abs(coord1.y - coord2.y);
        double z = cycle_abs(coord1.z - coord2.z);

        return sqrt(x*x + y*y + z*z);
    }

    // int getNumAtomInBasis(Vector atom) {
    //     atom.modf();
    //     int index = 0;
    //     for (const Vector& basisAtom : fccBasis) {
    //         if (atom.equals(basisAtom)) {
    //             break;
    //         }
    //         index++;
    //     }
    //     return index;
    // }

    // int getIndex(Vector& atom) {
    //     int n = getNumAtomInBasis(atom);
    //     int i = static_cast<int>(atom.x);
    //     int j = static_cast<int>(atom.y);
    //     int k = static_cast<int>(atom.z);

    //     return i + j*N + k*N*N + n*N*N*N;
    // }
};

struct Energy {
    double xi = 1.178;
    double r0 = 4.085 / sqrt(2);
    double q = 3.139;
    double A1 = 0, A0 = 0.1028, p = 10.928;

    double computeBindingEnergy(double r) {
        double tmp = r/r0 - 1;
        return exp(-2 * q * tmp);
    }

    double computeRepulsionEnergy(double r) {
        double tmp = r/r0 - 1;
        return (A1 * (r - r0) + A0) * exp(-1 * p * tmp);
    }

    double computeEnergy(Lattice& lt) {
        double Eb, Er, E = 0;
        int i = 0;
        for (Atom atom : lt.atoms) {
            //cout<< i <<": ";
            i++;
            int j = 0;
            Eb = Er = 0;
            for (Atom nb : lt.atoms) {
                double relative_radius = lt.dist(atom, nb);
                if (relative_radius > lt.cutoff || relative_radius == 0) {
                //     if (i == 1) {
                //     nb.coord.print();
                //     cout << endl;
                // }
                    cout << relative_radius << endl;
                    continue;
                }
                    
                double r = lt.a * relative_radius;
                // cout << j << ", ";
                // j++;
                if (i == 1) {
                    nb.coord.print();
                    cout<< ", " << endl;
                }
                Eb += computeBindingEnergy(r);
                Er += computeRepulsionEnergy(r);
            }
            E += Er - xi*sqrt(Eb);
            cout << endl << "E = " << E << endl;
        }
        return E;
    }
};

int main()
{
    int N = 3;
    int K = 4*N*N*N;
    Lattice lattice = Lattice(N, 4.085);
    Energy energy = Energy();
    cout << energy.computeEnergy(lattice)/K << endl;


    // ofstream outputFile("data.txt");
    // if (outputFile.is_open()) {
    //     double a = 2;
    //     cout << "равномерно" << endl;
    //     for (int i = 0; i < 100; i++) {
    //         a +=  0.01 * i;
    //         lattice.a = a;
    //         double y = energy.computeEnergy(lattice);
    //         outputFile << a  << " " << y << "\n";
    //     }
    //     outputFile.close();
    // } else {
    //     cerr << "Не удалось открыть файл для записи!" << endl;
    // }

    return 0;
}