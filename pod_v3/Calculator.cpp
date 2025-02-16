#include "definitions.h"

Calculator::Calculator(Lattice& lattice_init, ParametersGroup& pg_init, double delta_init) 
: lattice(lattice_init), params(pg_init), delta(delta_init) {}


//////////////////////////////////////////// Рассчет энергии ///////////////////////////////////////////////

Interaction Calculator::get_interaction(const Atom& atom1, const Atom& atom2) const {
    Metal metal1 = atom1.type;
    Metal metal2 = atom2.type;
    if (metal1 == metal2) {
        return (metal1 == Metal::ARGENTUM) ? Interaction::BASE : Interaction::EXTRA;
    } 
    return Interaction::SOLUTION;
}

Vector Calculator::get_mod_dist_vector(const Atom& atom1, const Atom& atom2, bool vacuum) const {
    Vector dist_vector = atom1.coord - atom2.coord;
    return dist_vector.mod(this->lattice.N, vacuum); 
}

double Calculator::compute_binding_energy(double r, Interaction interaction) const {
    Parameters p;
    switch (interaction) {
        case Interaction::BASE:
            p = this->params.paramAA;
            break;
        case Interaction::SOLUTION:
            p = this->params.paramAB;
            break;
        case Interaction::EXTRA:
            p = this->params.paramBB;
            break;
    }
    double tmp = r/p.r0 - 1;
    return p.xi* p.xi * exp(-2 * p.q * tmp);
}

double Calculator::compute_repulsion_energy(double r, Interaction interaction) const {
    Parameters p;
    switch (interaction) {
        case Interaction::BASE:
            p = this->params.paramAA;
            break;
        case Interaction::SOLUTION:
            p = this->params.paramAB;
            break;
        case Interaction::EXTRA:
            p = this->params.paramBB;
            break;
    }
    double tmp = r/p.r0 - 1;
    return (p.A1 * (r - p.r0) + p.A0) * exp(-1 * p.pr * tmp);
}

double Calculator::compute_energy(bool vacuum) const {
    double E = 0.0;

    // #pragma omp parallel for reduction(+:E) // Создаем параллельный цикл с редукцией
    for (size_t i = 0; i < this->lattice.atoms.size(); ++i) {
        const Atom& atom1 = this->lattice.atoms[i];
        double Eb = 0.0, Er = 0.0;

        for (const Atom& atom2 : this->lattice.atoms) {
            if (atom1 == atom2) continue;

            Interaction interaction = get_interaction(atom1, atom2);
            Vector mod_dist_vector = get_mod_dist_vector(atom1, atom2, vacuum);

            if (mod_dist_vector.length() > this->lattice.cutoff) continue;

            Vector norm_dist_vector = mod_dist_vector * this->lattice.a;
            double r = norm_dist_vector.length();

            Eb += compute_binding_energy(r, interaction);
            Er += compute_repulsion_energy(r, interaction);
        }

        // #pragma omp critical // Критическая секция для обновления общей переменной
        // {
            E += Er - sqrt(Eb);
        // }
    }

    return E;
}

double Calculator::compute_energy(vector<vector<double>> &deformation_matrix, bool vacuum) const {
    Vector a_prev = Vector(lattice.a); // точно копирует?
    lattice.a = lattice.a * deformation_matrix;
    double E = compute_energy(vacuum);
    lattice.a = Vector(a_prev);
    return E;
}

double Calculator::compute_cohesion_energy() const {
    return compute_energy()/lattice.K;
}


//////////////////////////////////////////// Рассчет констант упругости ///////////////////////////////////////////////

double Calculator::derivative(vector<vector<double>> &deformation_matrix_minus, vector<vector<double>> &deformation_matrix_plus) {
    double E1 = compute_energy(deformation_matrix_minus);
    double E2 = compute_energy();
    double E3 = compute_energy(deformation_matrix_plus);

    double second_derivative = (E1 - 2*E2 + E3)/(this->delta*this->delta);
    return second_derivative;
}

double Calculator::derivative_B() {
    vector<vector<double>> deformation_matrix_minus = {
        {1 - delta, 0.0, 0.0},
        {0.0, 1 - delta, 0.0},
        {0.0, 0.0, 1 - delta}
    };

    vector<vector<double>> deformation_matrix_plus = {
        {1 + delta, 0.0, 0.0},
        {0.0, 1 + delta, 0.0},
        {0.0, 0.0, 1 + delta}
    };

    return derivative(deformation_matrix_minus, deformation_matrix_plus);
}

double Calculator::derivative_C_11() {
    vector<vector<double>> deformation_matrix_minus = {
        {1 - delta, 0.0, 0.0},
        {0.0, 1 - delta, 0.0},
        {0.0, 0.0, 1.0}
    };

    vector<vector<double>> deformation_matrix_plus = {
        {1 + delta, 0.0, 0.0},
        {0.0, 1 + delta, 0.0},
        {0.0, 0.0, 1.0}
    };

    return derivative(deformation_matrix_minus, deformation_matrix_plus);
}

double Calculator::derivative_C_12() {
    vector<vector<double>> deformation_matrix_minus = {
        {1 - delta, 0.0, 0.0},
        {0.0, 1 + delta, 0.0},
        {0.0, 0.0, 1.0}
    };

    vector<vector<double>> deformation_matrix_plus = {
        {1 + delta, 0.0, 0.0},
        {0.0, 1 - delta, 0.0},
        {0.0, 0.0, 1.0}
    };

    return derivative(deformation_matrix_minus, deformation_matrix_plus);
}

double Calculator::derivative_C_44() {
    vector<vector<double>> deformation_matrix_minus = {
        {1.0, 0-delta, 0.0},
        {0-delta, 1.0, 0.0},
        {0.0, 0.0, 1/(1-delta*delta)}
    };

    vector<vector<double>> deformation_matrix_plus = {
        {1.0, delta, 0.0},
        {delta, 1.0, 0.0},
        {0.0, 0.0, 1/(1-delta*delta)}
    };

    return derivative(deformation_matrix_minus, deformation_matrix_plus);
}

ElasticityValues Calculator::get_elasticity() {
    double deriv_B = derivative_B();
    double deriv_C_11 = derivative_C_11();
    double deriv_C_12 = derivative_C_12();
    double deriv_C_44 = derivative_C_44();

    double B_value = (2*deriv_B)/(9*this->lattice.V);
    double C_11 = (deriv_C_11 + deriv_C_12)/(this->lattice.V * 2);
    double C_12 = (deriv_C_11 - deriv_C_12)/(this->lattice.V * 2);
    double C_44 = deriv_C_44/(this->lattice.V * 4);

    double d = 0.624*2;

    return {B_value/d, C_11/d, C_12/d, C_44/d};
}

void Calculator::elasticity_print() {
    ElasticityValues elasticity_values = get_elasticity();

    double C_11 = elasticity_values.C_11;
    double C_12 = elasticity_values.C_12;
    double B_value = elasticity_values.B;
    double C_44 = elasticity_values.C_44;


    cout << "B = " << B_value << endl;
    cout << "C11 = " << C_11 << endl;
    cout << "C12 = " << C_12 << endl;
    cout << "C44 = " << C_44 << endl;
}

//////////////////////////////////////////// Рассчет энергии растворимости ///////////////////////////////////////////////

double Calculator::compute_solubility_energy() {
    double E_B = compute_energy();
    double E_B_coh = E_B/this->lattice.K;
    double E_A_coh = -5.31; // когезионная энергия из таблицы для ванадия

    this->lattice.atoms[0].type = Metal::VANADIUM;
    double E_AB = compute_energy();
    this->lattice.atoms[0].type = Metal::ARGENTUM;

    double E_sol = E_AB - E_B - E_A_coh + E_B_coh;

    return E_sol;
}

//////////////////////////////////////////// Рассчет энергии связи димера в поверхностном слое ///////////////////////////////////////////////

double Calculator::compute_energy_of_dimer_in_surface() {
    double E_surf = compute_energy(true);

    int idx1 = this->lattice.get_index(Vector(1,1,0));
    int idx2 = this->lattice.get_index(Vector(1.5,1.5,0));

    this->lattice.atoms[idx1].type = Metal::VANADIUM;
    double E_adatom_surf = compute_energy(true);

    this->lattice.atoms[idx2].type = Metal::VANADIUM;
    double E_dim_surf = compute_energy(true);

    this->lattice.atoms[idx2].type = Metal::ARGENTUM;
    this->lattice.atoms[idx1].type = Metal::ARGENTUM;

    double E_dim = (E_dim_surf - E_surf) - 2*(E_adatom_surf - E_surf);

    return E_dim;

}

//////////////////////////////////////////// Рассчет энергии связи димера на поверхности ///////////////////////////////////////////////

double Calculator::compute_energy_of_dimer_on_surface() {
    double E_surf = compute_energy(true);

    Atom atom1 = Atom(Vector(1.5, 2, -0.5), Metal::VANADIUM);
    Atom atom2 = Atom(Vector(2, 2.5, -0.5), Metal::VANADIUM);

    this->lattice.atoms.push_back(atom1);
    double E_adatom_surf = compute_energy(true);

    this->lattice.atoms.push_back(atom2);
    double E_dim_surf = compute_energy(true);

    double E = compute_energy(true);

    this->lattice.atoms.pop_back();
    this->lattice.atoms.pop_back();

    double E_dim = (E_dim_surf - E_surf) - 2*(E_adatom_surf - E_surf);

    return E_dim;
}
