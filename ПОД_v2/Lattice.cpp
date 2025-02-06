#include "classes.h"
#include <iostream>
#include <cmath>

 // Стандартный конструктор
Lattice::Lattice(int N_init, double a_init, double cutoff_init) {
    N = N_init;
    K = 4*N*N*N;
    a = Vector(a_init, a_init, a_init);
    V = a_init * a_init * a_init * N * N * N;
    cutoff = cutoff_init;
    neighbours = initNeighbours();
    atoms = initLattice(N, Metal::ARGENTUM);
}

    // Конструктор для деформированной матрицы (есть копирование атомов и соседей) + (не регулируем cutoff)
Lattice::Lattice(const Lattice& original, vector<vector<double>> deltaMatrix) {
    N = original.N;
    K = original.K;
    a = original.a.multiplyByMatrix(deltaMatrix);
    V = a.x * a.y * a.z * N * N * N;
    cutoff = original.cutoff;
    fccBasis = original.fccBasis;
    neighbours = initNeighbours();
    atoms = initLattice(N, Metal::ARGENTUM);
}


bool Lattice::isAtomInFCC(Vector atom) {
    atom.modf();
    for (const Vector& basisAtom : fccBasis) {
        if (atom.equals(basisAtom)) {
            return true;
        }
    }

    return false;
}

vector<Vector> Lattice::initNeighbours() {
    vector<Vector> neighbours;
    double coordCutoff = static_cast<int>(cutoff / 0.5) * 0.5;
    Vector center = Vector(0, 0, 0);
    Vector nb = Vector(center);
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

vector<Atom> Lattice::initLattice(int N, Metal metal) {
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


int Lattice::getIndexAtom(Vector v) const{
    int idx = 0;
    for (const Atom atom : atoms) {
        Vector atom_coord = atom.coord;
        if(atom_coord==v) {
            return idx;
        }
        idx++;
    }
    return -1;
}


