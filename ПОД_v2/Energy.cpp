#include "classes.h"
#include <iostream>
#include <omp.h>



Energy::Energy() {
    pg  = ParametersGroup();
}

Energy::Energy(ParametersGroup &paramGroup) {
    pg  = ParametersGroup(paramGroup);
}

double Energy::computeBindingEnergy(double r, Structure s) const {
    Parameters p;
    switch (s) {
        case Structure::BASE:
            p = pg.paramAA;
            break;
        case Structure::SOLUTION:
            p = pg.paramAB;
            break;
        case Structure::EXTRA:
            p = pg.paramBB;
            break;
    }
    double tmp = r/p.r0 - 1;
    return p.xi* p.xi * exp(-2 * p.q * tmp);
}

double Energy::computeRepulsionEnergy(double r, Structure s) const {
    Parameters p;
    switch (s) {
        case Structure::BASE:
            p = pg.paramAA;
            break;
        case Structure::SOLUTION:
            p = pg.paramAB;
            break;
        case Structure::EXTRA:
            p = pg.paramBB;
            break;
    }
    double tmp = r/p.r0 - 1;
    return (p.A1 * (r - p.r0) + p.A0) * exp(-1 * p.pr * tmp);
}

double Energy::computeEnergy(const Lattice& lt, bool vacuum) const {
    double lower_limit = 0, upper_limit = lt.a.length() * lt.N;
    double E = 0;

    #pragma omp parallel for reduction(+:E)
    for (size_t i = 0; i < lt.atoms.size(); ++i) {
        const Atom& atom = lt.atoms[i];
        double Eb = 0, Er = 0;
        Metal metal1 = atom.type;

        for (const Vector& nb_coord : lt.neighbours) {
            if (vacuum && (nb_coord.z + atom.coord.z < lower_limit || nb_coord.z + atom.coord.z > upper_limit)) {
                continue;
            }

            Structure s;
            int nb_idx = lt.getIndexAtom(nb_coord);
            if (nb_idx == -1) {
                s = (metal1 == Metal::ARGENTUM) ? Structure::BASE : Structure::EXTRA;
            } else {
                const Atom& nb = lt.atoms[nb_idx];
                Metal metal2 = nb.type;
                if (metal1 == metal2) {
                    s = (metal1 == Metal::ARGENTUM) ? Structure::BASE : Structure::EXTRA;
                } else {
                    s = Structure::SOLUTION;
                }
            }

            double r = nb_coord.multiplyByVector(lt.a).length();
            Eb += computeBindingEnergy(r, s);
            Er += computeRepulsionEnergy(r, s);
        }
        #pragma omp critical
        E += Er - sqrt(Eb); 
    }
    return E;
}

double Energy::computeEnergyOne(double r, Structure s) const {
    double E = 0;

    double Eb = computeBindingEnergy(r, s);
    double Er = computeRepulsionEnergy(r, s);
    E += Er - sqrt(Eb); 
    return E;
}

double Energy::computeCohesionEnergy(const Lattice& lt) const {
    return computeEnergy(lt)/lt.K;
}

double Energy::computeSolubilityEnergy(Lattice& lt) {
    double E_B = computeEnergy(lt);
    double E_B_coh = E_B/lt.K;
    double E_A_coh = -4.28;

    lt.atoms[0].type = Metal::FERRUM;
    double E_AB = computeEnergy(lt);
    lt.atoms[0].type = Metal::ARGENTUM;

    double E_sol = E_AB - E_B - E_A_coh + E_B_coh;

    return E_sol;
}



double Energy::computeEnergyOfDimerInSurface(Lattice& lt) {
    double E_surf = computeEnergy(lt, true);
    int n = 2;


    int nb_idx1 = lt.getIndexAtom(Vector(1,1,0));
    lt.atoms[nb_idx1].type = Metal::FERRUM;
    double E_adatom_surf = computeEnergy(lt, true);

    int nb_idx2 = lt.getIndexAtom(Vector(1.5,1.5,0));
    lt.atoms[nb_idx2].type = Metal::FERRUM;
    double E_dim_surf = computeEnergy(lt, true);

    lt.atoms[nb_idx2].type = Metal::ARGENTUM;
    lt.atoms[nb_idx1].type = Metal::ARGENTUM;

    double E_dim = (E_dim_surf - E_surf) - n*(E_adatom_surf - E_surf);

    return E_dim;
}

double Energy::computeEnergyOfDimerOnSurface(Lattice &lattice) {
    Lattice lt_more = Lattice(4, lattice.a.length());
    double E_surf = computeEnergy(lt_more, true);
    int n = 2;

    double Eb1 = 0, Er1 = 0;
    double Eb2 = 0, Er2 = 0;
    Atom atom1 = Atom(Vector(1.5, 2, -0.5), Metal::FERRUM);
    Atom atom2 = Atom(Vector(2, 2.5, -0.5), Metal::FERRUM);
    Structure s = Structure::SOLUTION;

    // Параллелизация цикла
    #pragma omp parallel for reduction(+:Eb1, Er1, Eb2, Er2)
    for (size_t i = 0; i < lt_more.atoms.size(); ++i) {
        Atom atom = lt_more.atoms[i];
        
        Vector r1_coord = atom.coord - atom1.coord;
        Vector r2_coord = atom.coord - atom2.coord;

        if (r1_coord.length() < lt_more.cutoff) {
            int k1 = (atom.coord.x >= 3 || atom.coord.y < 1) ? 1 : 2;
            double r1 = r1_coord.multiplyByVector(lt_more.a).length();
            Eb1 += k1 * computeBindingEnergy(r1, s);
            Er1 += k1 * computeRepulsionEnergy(r1, s);
        }

        if (r2_coord.length() < lt_more.cutoff) {
            int k2 = (atom.coord.x >= 3 || atom.coord.y < 1) ? 1 : 2;
            double r2 = r2_coord.multiplyByVector(lt_more.a).length();
            Eb2 += k2 * computeBindingEnergy(r2, s);
            Er2 += k2 * computeRepulsionEnergy(r2, s);
        }
    }

    double E_adatom_surf = Er1 - sqrt(Eb1) + E_surf;

    double r3 = (atom2.coord - atom1.coord).multiplyByVector(lt_more.a).length();
    Structure s3 = Structure::EXTRA;
    double Eb3 = computeBindingEnergy(r3, s3);
    double Er3 = computeRepulsionEnergy(r3, s3);
    double E_dim_surf = (Er1 + Er3 - sqrt(Eb1 + Eb3)) + (Er2 + Er3 - sqrt(Eb2 + Eb3)) + E_surf;

    double E_dim = (E_dim_surf - E_surf) - n * (E_adatom_surf - E_surf);

    return E_dim;
}