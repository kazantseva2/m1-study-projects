#include "definitions.h"


Lattice::Lattice(int N_init, double a_init, double cutoff_init) {
    N = N_init;
    K = 4*N*N*N;
    a = Vector(a_init, a_init, a_init);
    V = a_init * a_init * a_init * N * N * N;
    cutoff = cutoff_init;
    atoms = init_lattice(N, Metal::ARGENTUM);
}

vector<Atom> Lattice::init_lattice(int N, Metal metal) {
    vector<Atom> atoms;
    for (Vector basisAtom : fcc_basis) {
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

int Lattice::get_index(const Vector& coord) const {
    for (size_t i = 0; i < this->atoms.size(); i++) {
        if (atoms[i].coord == coord) {
            return i;
        }
    }
    return -1; 
}