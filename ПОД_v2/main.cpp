#include <iostream>
#include <chrono>
#include "classes.h"
#include <fstream>

// Ctrl+Shift+B
// g++ -fopenmp main.cpp Vector.cpp Lattice.cpp Energy.cpp Elasticity.cpp Atom.cpp Optimize.cpp Parameters.cpp -o main
// ./main

int main() {
    int N = 1;
    double a = 4.085;
    Lattice lattice = Lattice(N, a);
    Energy energy = Energy();
    cout << energy.computeEnergy(lattice) << endl;
    //cout << energy.computeCohesionEnergy(lattice) << endl;

    // Parameters paramAA = Parameters();
    // paramAA.A0 = 0.2061; paramAA.A1 = 0.0; paramAA.pr = 10.229; paramAA.q = 4.036;
    // paramAA.xi = 1.790; paramAA.r0 = a/sqrt(2);

    // // paramAA.A0 = 0.0992012; paramAA.A1 = -0.0540063; paramAA.pr = 9.02197; paramAA.q = 3.00369;
    // // paramAA.xi = 1.1683; paramAA.r0 = 2.84484;

    // Parameters paramAB = Parameters();
    // // paramAB.A0 = 0.0873938; paramAB.A1 = -0.0361787; paramAB.pr = 12.3725; paramAB.q = 2.11971;
    // // paramAB.xi = 1.31401; paramAB.r0 = 2.6801;

    // Parameters paramBB = Parameters();
    // // paramBB.A0 = 0.102725; paramBB.A1 = -0.0542089; paramBB.pr = 11.1873; paramBB.q = 3.38134;
    // // paramBB.xi = 1.33753; paramBB.r0 = 2.65177;

    // ParametersGroup pg = ParametersGroup();
    // pg.paramAA = paramAA;
    // pg.paramAB = paramAB;
    // pg.paramBB = paramBB;

     // A-A : A0: 0.0992012 A1: -0.0540063 p0: 9.02197 q0: 3.00369 qsi: 1.1683 r0: 2.84484
    // A-B : A0: 0.0873938 A1: -0.0361787 p0: 12.3725 q0: 2.11971 qsi: 1.31401 r0: 2.6801
    // B-B : A0: 0.102725 A1: -0.0542089 p0: 11.1873 q0: 3.38134 qsi: 1.33753 r0: 2.65177
    // params : a=4.008, e_c=-2.99045, B=1.03718, C_11=1.24766, C_12=0.931918, C_44=0.481209,
    // e_sol=0.552351, e_in_dim=0.184095, e_on_dim=-0.321981

    // Energy energy = Energy(pg);
    // double E_coh = energy.computeCohesionEnergy(lattice);
    // Elasticity elasticity = Elasticity();
    // ElasticityValues elVls = elasticity.getElasticity(energy, lattice, 0.01);

    // double E_sol = energy.computeSolubilityEnergy(lattice);

    // double E_in_dim = energy.computeEnergyOfDimerInSurface(lattice);
    // double E_on_dim = energy.computeEnergyOfDimerOnSurface(lattice);


    // cout << "Значения оптимизации:" << endl;
    // cout << "  Параметр решетки a: " << a << endl;
    // cout << "  Когезионная энергия E_coh: " << E_coh  << endl;
    // cout << "  Упругость:" << endl;
    // cout << "    B: " << elVls.B << endl;
    // cout << "    C11: " << elVls.C11 << endl;
    // cout << "    C12: " << elVls.C12 << endl;
    // cout << "    C44: " << elVls.C44 << endl;


    // Parameters paramAA = Parameters();
    // paramAA.A0 = 0.0992012; paramAA.A1 = -0.0540063; paramAA.pr = 9.02197; paramAA.q = 3.00369;
    // paramAA.xi = 1.1683; paramAA.r0 = 2.84484;

    // Parameters paramAB = Parameters();
    // paramAB.A0 = 0.0873938; paramAB.A1 = -0.0361787; paramAB.pr = 12.3725; paramAB.q = 2.11971;
    // paramAB.xi = 1.31401; paramAB.r0 = 2.6801;

    // Parameters paramBB = Parameters();
    // paramBB.A0 = 0.102725; paramBB.A1 = -0.0542089; paramBB.pr = 11.1873; paramBB.q = 3.38134;
    // paramBB.xi = 1.33753; paramBB.r0 = 2.65177;

    // ParametersGroup pg = ParametersGroup();
    // pg.paramAA = paramAA;
    // pg.paramAB = paramAB;
    // pg.paramBB = paramBB;

    //  // A-A : A0: 0.0992012 A1: -0.0540063 p0: 9.02197 q0: 3.00369 qsi: 1.1683 r0: 2.84484
    // // A-B : A0: 0.0873938 A1: -0.0361787 p0: 12.3725 q0: 2.11971 qsi: 1.31401 r0: 2.6801
    // // B-B : A0: 0.102725 A1: -0.0542089 p0: 11.1873 q0: 3.38134 qsi: 1.33753 r0: 2.65177
    // // params : a=4.008, e_c=-2.99045, B=1.03718, C_11=1.24766, C_12=0.931918, C_44=0.481209,
    // // e_sol=0.552351, e_in_dim=0.184095, e_on_dim=-0.321981

    // Energy energy = Energy(pg);
    // double E_coh = energy.computeCohesionEnergy(lattice);
    // Elasticity elasticity = Elasticity();
    // ElasticityValues elVls = elasticity.getElasticity(energy, lattice, 0.01);

    // double E_sol = energy.computeSolubilityEnergy(lattice);

    // double E_in_dim = energy.computeEnergyOfDimerInSurface(lattice);
    // double E_on_dim = energy.computeEnergyOfDimerOnSurface(lattice);


    // cout << "Значения оптимизации:" << endl;
    // cout << "  Параметр решетки a: " << a << endl;
    // cout << "  Когезионная энергия E_coh: " << E_coh  << endl;
    // cout << "  Упругость:" << endl;
    // cout << "    B: " << elVls.B << endl;
    // cout << "    C11: " << elVls.C11 << endl;
    // cout << "    C12: " << elVls.C12 << endl;
    // cout << "    C44: " << elVls.C44 << endl;
    // cout << "  Энергия растворимости E_sol: " << E_sol  << endl;
    // cout << "  Энергия связи димера в поверхностном слое E_in_dim: " << E_in_dim << endl;
    // cout << "  Энергия связи димера на поверхностном слое E_on_dim: " << E_on_dim << endl;

    // // Засекаем время
    // auto start = std::chrono::high_resolution_clock::now();

    // Optimize opt = Optimize();
    // opt.optimize(lattice, 6, Structure::BASE, 0.001);

    // Energy energy1 = Energy(opt.pg_opt);
    // cout << "Оптимальные параметры AA:" << endl;
    // opt.pg_opt.paramAA.print();
   
    // double E_coh = energy1.computeCohesionEnergy(lattice);
    // Elasticity elasticity = Elasticity();
    // ElasticityValues elVls = elasticity.getElasticity(energy1, lattice, 0.01);

    // cout << "Значения оптимизации:" << endl;
    // cout << "  Параметр решетки a: " << a << " (opt: " << opt.a_opt << ")" << endl;
    // cout << "  Когезионная энергия E_coh: " << E_coh << " (opt: " << opt.E_coh_opt << ")" << endl;
    // cout << "  Упругость:" << endl;
    // cout << "    B: " << elVls.B << " (opt: " << opt.B_opt << ")" << endl;
    // cout << "    C11: " << elVls.C11 << " (opt: " << opt.C11_opt << ")" << endl;
    // cout << "    C12: " << elVls.C12 << " (opt: " << opt.C12_opt << ")" << endl;
    // cout << "    C44: " << elVls.C44 << " (opt: " << opt.C44_opt << ")" << endl;

    
    // opt.optimize(lattice, 6, Structure::SOLUTION, 0.0001);
    // Energy energy2 = Energy(opt.pg_opt);

    // cout << "Оптимальные параметры AB:" << endl;
    // opt.pg_opt.paramAB.print();

    // cout << "Значения оптимизации:" << endl;
    // double E_sol = energy2.computeSolubilityEnergy(lattice);
    // cout << "  Энергия растворимости E_sol: " << E_sol << " (opt: " << opt.E_sol_opt << ")" << endl;



    // opt.optimize(lattice, 6, Structure::EXTRA, 0.001);
    // Energy energy3 = Energy(opt.pg_opt);

    // cout << "Оптимальные параметры BB:" << endl;
    // opt.pg_opt.paramBB.print();

    // cout << "Значения оптимизации:" << endl;

    // double E_in_dim = energy3.computeEnergyOfDimerInSurface(lattice);
    // double E_on_dim = energy3.computeEnergyOfDimerOnSurface(lattice);
    // cout << "  Энергия связи димера в поверхностном слое E_in_dim: " << E_in_dim << " (opt: " << opt.E_in_dim_opt << ")" << endl;
    // cout << "  Энергия связи димера на поверхностном слое E_on_dim: " << E_on_dim << " (opt: " << opt.E_on_dim_opt << ")" << endl;


    // // Останавливаем таймер
    // auto end = std::chrono::high_resolution_clock::now();

    // // Вычисляем время в миллисекундах
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "Время работы: " << 1.0*duration.count()/1000 << " s" << std::endl;
   
   

    // ofstream outputFile("data.txt");
    // if (outputFile.is_open()) {
    //     cout << "равномерно" << endl;
    //     for (int i = 1; i < 10; i+=0.1) {
    //         double y = energy1.computeEnergyOne(i, Structure::BASE);
    //         outputFile << i  << " " << y << "\n";
    //     }
    //     outputFile.close();
    // } else {
    //     cerr << "Не удалось открыть файл для записи!" << endl;
    // }
    return 0;
}