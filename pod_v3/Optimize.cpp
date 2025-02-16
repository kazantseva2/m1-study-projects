#include "definitions.h"

Optimize::Optimize(Calculator& calculator_init) : calculator(calculator_init){
}

double Optimize::mse(double f_val, double f_opt, double w) {
    double error = f_val - f_opt;
    double res = w*error*error/abs(f_opt);
    return res;
}

double Optimize::error_AA(Parameters &p) {
    this->calculator.params.paramAA = p;
    double a = 4.085;
    double E_coh = this->calculator.compute_cohesion_energy();
    ElasticityValues elasticity_values = this->calculator.get_elasticity();
    
    double res = 0;
    double w = 1.0/6;
    res += mse(a, a_opt, w);
    res += mse(E_coh, E_coh_opt, w);
    res += mse(elasticity_values.B, this->B_opt, w);
    res += mse(elasticity_values.C_11, this->C_11_opt, w);
    res += mse(elasticity_values.C_12, this->C_12_opt, w);
    res += mse(elasticity_values.C_44, this->C_44_opt, w);

    p.error_value = res;
    return res;
}

double Optimize::error_AB(Parameters &p) {
    this->calculator.params.paramAB = p;
    double E_sol = this->calculator.compute_solubility_energy();
    
    double res = mse(E_sol, this->E_sol_opt, 1.0);

    p.error_value = res;

    return res;
}

double Optimize::error_BB(Parameters &p) {
    this->calculator.params.paramBB = p;
    
    double E_in_dim = this->calculator.compute_energy_of_dimer_in_surface();
    double E_on_dim = this->calculator.compute_energy_of_dimer_on_surface();
    
    double res = 0;
    double w = 1.0/2;
    res += mse(E_in_dim, this->E_in_dim_opt, w);
    res += mse(E_on_dim, this->E_on_dim_opt, w);

    p.error_value = res;

    return res;
}

void Optimize::error_fun(Parameters &p, Interaction interact) {
    switch (interact) {
        case Interaction::BASE:
            error_AA(p);
            break;
        case Interaction::SOLUTION:
            error_AB(p);
            break;
        case Interaction::EXTRA:
            error_BB(p);
            break;
    }
}

void Optimize::initialize_simplex(Interaction interact, int n) {
    simplex.clear();

    for (int i = 0; i < n+1; i++) {
        Parameters param;
        
        switch (interact) {
            case Interaction::BASE:
                param.rand_init();
                param.error_value = error_AA(param);
                break;
            case Interaction::SOLUTION:
                param.rand_init();
                param.error_value = error_AB(param);
                break;
            case Interaction::EXTRA:
                param.rand_init();
                param.error_value = error_BB(param);
                break;
        }
        simplex.push_back(param);
    }
}

void Optimize::sort_simplex_by_value() {
    sort(simplex.begin(), simplex.end(), [](const Parameters& a, const Parameters& b) {
        return a.error_value > b.error_value; // Сортировка по убыванию
    });
}

Parameters Optimize::get_center() {
    Parameters x_c; 
    int n = simplex.size();

    x_c.xi = 0.0;
    x_c.r0 = 0.0;
    x_c.q = 0.0;
    x_c.A1 = 0.0;
    x_c.A0 = 0.0;
    x_c.pr = 0.0;

    for (int i = 1; i < n; i++) {
        x_c += simplex[i];
    }

    double scale = 1.0 / (n - 1);
    x_c *= scale;

    return x_c;
}

void Optimize::optimize(Interaction interact, int n, double eps) {

    initialize_simplex(interact, n);

    for (int i = 0; i < 100; i++) {
        sort_simplex_by_value();
        bool compression = false;
        Parameters x_c = get_center();


        // отражение
        Parameters x_r =  x_c * (1 + alpha)  -  simplex[0] * alpha;
        error_fun(x_r, interact);

        if (x_r.error_value < simplex[n].error_value) { 
            Parameters x_e =  x_c * (1 - gamma)  +  x_r * gamma;
            error_fun(x_e, interact);
            if (x_e.error_value < x_r.error_value) {
                simplex[0] = x_e;
            } else {
                simplex[0] = x_r;
            }
        } else if (x_r.error_value < simplex[1].error_value) {
            simplex[0] = x_r;
        } else if (x_r.error_value < simplex[0].error_value) {
            compression = true;
            simplex[0] = x_r;
        } else {
            compression = true;
        }

        // сжатие
        if (compression) {
            Parameters x_s = x_c * (1 - beta)  +  simplex[0] * beta;
            error_fun(x_s, interact);
            if (x_s.error_value < simplex[0].error_value) {
                simplex[0] = x_s;
            } else {
                for (int i = 0; i < n; i++) {
                    simplex[i] = simplex[n] + (simplex[i] - simplex[n]) * sigma;
                    error_fun(simplex[i], interact);
                }
            }
            compression = false;
        }

        if (simplex[n].error_value < eps) {
            break; // Алгоритм сошелся
        }
    }
    
    sort_simplex_by_value();
    switch (interact) {
        case Interaction::BASE:
            this->calculator.params.paramAA = simplex[n];
            break;
        case Interaction::SOLUTION:
            this->calculator.params.paramAB = simplex[n];
            break;
        case Interaction::EXTRA:
            this->calculator.params.paramBB = simplex[n];
            break;
    }
}

void Optimize::print_parameters() const {
    cout << "Оптимальные параметры AA:" << endl;
    this->calculator.params.paramAA.print();
    cout << "Оптимальные параметры AB:" << endl;
    this->calculator.params.paramAB.print();
    cout << "Оптимальные параметры BB:" << endl;
    this->calculator.params.paramBB.print();
}

void Optimize::print_result() const {
    print_parameters();

    double a = 4.085;
    double E_coh = calculator.compute_cohesion_energy();
    ElasticityValues elasticity_values = calculator.get_elasticity();
    double E_sol = calculator.compute_solubility_energy();
    double E_in_dim = calculator.compute_energy_of_dimer_in_surface();
    double E_on_dim = calculator.compute_energy_of_dimer_on_surface();

    cout << "Значения оптимизации:" << endl;
    cout << "  Параметр решетки a: " << a << " (opt: " << this->a_opt << ")" << endl;
    cout << "  Когезионная энергия E_coh: " << E_coh << " (opt: " << this->E_coh_opt << ")" << endl;
    cout << "  Упругость:" << endl;
    cout << "    B: " << elasticity_values.B << " (opt: " << this->B_opt << ")" << endl;
    cout << "    C11: " << elasticity_values.C_11 << " (opt: " << this->C_11_opt << ")" << endl;
    cout << "    C12: " << elasticity_values.C_12 << " (opt: " << this->C_12_opt << ")" << endl;
    cout << "    C44: " << elasticity_values.C_44 << " (opt: " << this->C_44_opt << ")" << endl;
    cout << "  Энергия растворимости E_sol: " << E_sol << " (opt: " << this->E_sol_opt << ")" << endl;
    cout << "  Энергия связи димера в поверхностном слое E_in_dim: " << E_in_dim << " (opt: " << this->E_in_dim_opt << ")" << endl;
    cout << "  Энергия связи димера на поверхностном слое E_on_dim: " << E_on_dim << " (opt: " << this->E_on_dim_opt << ")" << endl;
}



