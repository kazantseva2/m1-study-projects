#include "definitions.h"
#include <chrono>

//  g++ -fopenmp main.cpp Vector.cpp Lattice.cpp Calculator.cpp Atom.cpp Optimize.cpp Parameters.cpp -o main
// ./main

//  g++ main.cpp Vector.cpp Lattice.cpp Calculator.cpp Atom.cpp Optimize.cpp Parameters.cpp -o main
// ./main

int main() {
    int N = 3;
    double a = 4.085;
    double delta = 0.01;
    Lattice lattice = Lattice(N, a);
    ParametersGroup param_group = ParametersGroup();
    Calculator calculator = Calculator(lattice, param_group, delta);


    cout << "Энергия: "<< calculator.compute_energy() << endl;
    cout << "Когезионная энергия: "<< calculator.compute_cohesion_energy() << endl;

    calculator.elasticity_print();

    cout << "Энергия растворимости: "<< calculator.compute_solubility_energy() << endl;
    cout << "Энергия связи димера в поверхностном слое: "<< calculator.compute_energy_of_dimer_in_surface() << endl;
    cout << "Энергия связи димера на поверхности: "<< calculator.compute_energy_of_dimer_on_surface() << endl;

    // Оптимизация
    Optimize optimization = Optimize(calculator);

    auto start = chrono::high_resolution_clock::now();
    optimization.optimize(Interaction::BASE, 6, 0.001);
    auto end = std::chrono::high_resolution_clock::now();
    cout << "BASE optimization time: " << chrono::duration<double>(end - start).count() << " seconds\n";

    start = std::chrono::high_resolution_clock::now();
    optimization.optimize(Interaction::SOLUTION, 2, 0.001);
    end = chrono::high_resolution_clock::now();
    cout << "SOLUTION optimization time: "  << chrono::duration<double>(end - start).count() << " seconds\n";

    start = chrono::high_resolution_clock::now();
    optimization.optimize(Interaction::EXTRA, 2, 0.001);
    end = chrono::high_resolution_clock::now();
    cout << "EXTRA optimization time: " << chrono::duration<double>(end - start).count() << " seconds\n";

    optimization.print_result();

    return 0;
}