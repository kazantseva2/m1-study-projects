#include "classes.h"
#include <iostream>


using namespace std;

double Optimize::mse(double f_val, double f_opt, double w) {
    double error = f_val - f_opt;
    double res = w*error*error/abs(f_opt);
    return res;
}

double Optimize::errorAA(Parameters &p, Lattice &lt) {
    ParametersGroup pg;
    pg.paramAA = p;
    Energy energy = Energy(pg);
    double a = 4.085;
    double E_coh = energy.computeCohesionEnergy(lt);

    Elasticity elasticity = Elasticity();
    ElasticityValues elVls = elasticity.getElasticity(energy, lt, 0.01);
    
    double res = 0;
    double w = 1.0/6;
    res += mse(a, a_opt, w);
    res += mse(E_coh, E_coh_opt, w);
    res += mse(elVls.B, B_opt, w);
    res += mse(elVls.C11, C11_opt, w);
    res += mse(elVls.C12, C12_opt, w);
    res += mse(elVls.C44, C44_opt, w);

    p.error = res;
    return res;
}

double Optimize::errorAB(Parameters &p, Lattice &lt) {
    ParametersGroup pg;
    pg.paramAA = pg_opt.paramAA;
    pg.paramAB = p;
    Energy energy = Energy(pg);
    double E_sol = energy.computeSolubilityEnergy(lt);
    
    double res = mse(E_sol, E_sol_opt, 1.0);

    p.error = res;

    return res;
}

double Optimize::errorBB(Parameters &p, Lattice &lt) {
    ParametersGroup pg;
    pg.paramAA = pg_opt.paramAA;
    pg.paramAB = pg_opt.paramAB;
    pg.paramBB = p;
    Energy energy = Energy(pg);
    
    double E_in_dim = energy.computeEnergyOfDimerInSurface(lt);
    double E_on_dim = energy.computeEnergyOfDimerOnSurface(lt);
    
    double res = 0;
    double w = 1.0/2;
    res += mse(E_in_dim, E_in_dim_opt, w);
    res += mse(E_on_dim, E_on_dim_opt, w);

    p.error = res;

    return res;
}

void Optimize::errorFun(Parameters &p, Lattice &lt, Structure st) {
    switch (st) {
            case Structure::BASE:
                errorAA(p, lt);
                break;
            case Structure::SOLUTION:
                errorAB(p, lt);
                break;
            case Structure::EXTRA:
                errorBB(p, lt);
                break;
        }
}

void Optimize::initializeSimplex(Lattice &lt, Structure st, int n) {
    simplex.clear();
    

    for (int i = 0; i < n+1; i++) {
        Parameters param;
        
        switch (st) {
            case Structure::BASE:
                param.randInit();
                param.error = errorAA(param, lt);
                break;
            case Structure::SOLUTION:
                param.randInit();
                param.error = errorAB(param, lt);
                break;
            case Structure::EXTRA:
                param.randInit();
                param.error = errorBB(param, lt);
                break;
        }
        simplex.push_back(param);
    }
}

void Optimize::sortSimplexByValue() {
    sort(simplex.begin(), simplex.end(), [](const Parameters& a, const Parameters& b) {
        return a.error > b.error; // Сортировка по убыванию
    });
}

Parameters Optimize::getCenter() {
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

    // F(x_c, lt);
    return x_c;
}

    
ParametersGroup Optimize::optimize(Lattice &lt, int n, Structure st, double eps) {

    initializeSimplex(lt, st);

    for (int i = 0; i < 100; i++) {
        sortSimplexByValue();
        bool compression = false;
        Parameters x_c = getCenter();


        // отражение
        Parameters x_r =  x_c * (1 + alpha)  -  simplex[0] * alpha;
        errorFun(x_r, lt, st);

        if (x_r.error < simplex[n].error) { 
            Parameters x_e =  x_c * (1 - gamma)  +  x_r * gamma;
            errorFun(x_e, lt, st);
            if (x_e.error < x_r.error) {
                simplex[0] = x_e;
            } else {
                simplex[0] = x_r;
            }
        } else if (x_r.error < simplex[1].error) {
            simplex[0] = x_r;
        } else if (x_r.error < simplex[0].error) {
            compression = true;
            simplex[0] = x_r;
        } else {
            compression = true;
        }

        // сжатие
        if (compression) {
            Parameters x_s = x_c * (1 - beta)  +  simplex[0] * beta;
            errorFun(x_s, lt, st);
            if (x_s.error < simplex[0].error) {
                simplex[0] = x_s;
            } else {
                for (int i = 0; i < n; i++) {
                    simplex[i] = simplex[n] + (simplex[i] - simplex[n]) * sigma;
                    errorFun(simplex[i], lt, st);
                }
            }
            compression = false;
        }

        // // проверка сходимости
        // double maxDifference = 0.0;
        // for (int i = 1; i < n; ++i) {
        //     maxDifference = max(maxDifference, std::abs(simplex[i].error - simplex[0].error));
        // }
        // cout << maxDifference << endl;
        // if (maxDifference < 0.00001) {
        //     break; // Алгоритм сошелся
        // }
        
        
        // cout << "error(" << i << ") = " << simplex[n].error << endl;
        if (simplex[n].error < eps) {
            break; // Алгоритм сошелся
        }
    }
    
    sortSimplexByValue();
    switch (st) {
        case Structure::BASE:
            pg_opt.paramAA = simplex[n];
            break;
        case Structure::SOLUTION:
            pg_opt.paramAB = simplex[n];
            break;
        case Structure::EXTRA:
            pg_opt.paramBB = simplex[n];
            break;
    }

    return pg_opt;
}

double Optimize::printResult(Lattice &lt) {
    Energy energy = Energy(pg_opt);
    cout << "Оптимальные параметры AA:" << endl;
    pg_opt.paramAA.print();
    cout << "Оптимальные параметры AB:" << endl;
    pg_opt.paramAB.print();
    cout << "Оптимальные параметры BB:" << endl;
    pg_opt.paramBB.print();

    double a = 4.085;
    
    double E_coh = energy.computeCohesionEnergy(lt);

    Elasticity elasticity = Elasticity();
    ElasticityValues elVls = elasticity.getElasticity(energy, lt, 0.01);

    double E_sol = energy.computeSolubilityEnergy(lt);
    double E_in_dim = energy.computeEnergyOfDimerInSurface(lt);
    double E_on_dim = energy.computeEnergyOfDimerOnSurface(lt);

    double res = errorAA(pg_opt.paramAA, lt) + errorAB(pg_opt.paramAA, lt) + errorBB(pg_opt.paramAA, lt);

    cout << "Значения оптимизации:" << endl;
    cout << "  Параметр решетки a: " << a << " (opt: " << a_opt << ")" << endl;
    cout << "  Когезионная энергия E_coh: " << E_coh << " (opt: " << E_coh_opt << ")" << endl;
    cout << "  Упругость:" << endl;
    cout << "    B: " << elVls.B << " (opt: " << B_opt << ")" << endl;
    cout << "    C11: " << elVls.C11 << " (opt: " << C11_opt << ")" << endl;
    cout << "    C12: " << elVls.C12 << " (opt: " << C12_opt << ")" << endl;
    cout << "    C44: " << elVls.C44 << " (opt: " << C44_opt << ")" << endl;
    cout << "  Энергия растворимости E_sol: " << E_sol << " (opt: " << E_sol_opt << ")" << endl;
    cout << "  Энергия связи димера в поверхностном слое E_in_dim: " << E_in_dim << " (opt: " << E_in_dim_opt << ")" << endl;
    cout << "  Энергия связи димера на поверхностном слое E_on_dim: " << E_on_dim << " (opt: " << E_on_dim_opt << ")" << endl;
    cout << "  Ошибка: " << res << endl;

    return res;
}



