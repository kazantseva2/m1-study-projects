#include "classes.h"
#include <iostream>


using namespace std;

void Parameters::print() const {
    cout << "Parameters:" << endl;

    cout << "  xi: ";
    cout << xi << " ";
    cout << endl;

    cout << "  r0: ";
    cout << r0 << " ";
    cout << endl;

    cout << "  q : ";
    cout << q << " ";
    cout << endl;

    cout << "  A1: ";
    cout << A1 << " ";
    cout << endl;

    cout << "  A0: ";
    cout << A0 << " ";
    cout << endl;

    cout << "  p : ";
    cout << pr << " ";
    cout << endl;
}

void Parameters::operator*=(double scalar) {
    xi *= scalar;
    r0 *= scalar;
    q *= scalar;
    A1 *= scalar;
    A0 *= scalar;
    pr *= scalar;
}

void Parameters::operator+=(Parameters other) {
    xi += other.xi;
    r0 += other.r0;
    q += other.q;
    A1 += other.A1;
    A0 += other.A0;
    pr += other.pr;
}

Parameters Parameters::operator*(double scalar) const {
    Parameters result;
    result.xi = xi * scalar;
    result.r0 = r0 * scalar;
    result.q = q * scalar;
    result.A1 = A1 * scalar;
    result.A0 = A0 * scalar;
    result.pr = pr * scalar;
    
    return result;
}

Parameters Parameters::operator-(Parameters other) const {
    Parameters result;
    result.xi = xi - other.xi;
    result.r0 = r0 - other.r0;
    result.q = q - other.q;
    result.A1 = A1 - other.A1;
    result.A0 = A0 - other.A0;
    result.pr = pr - other.pr;

    return result;
}

Parameters Parameters::operator+(Parameters other) const {
    Parameters result;
    result.xi = xi + other.xi;
    result.r0 = r0 + other.r0;
    result.q = q + other.q;
    result.A1 = A1 + other.A1;
    result.A0 = A0 + other.A0;
    result.pr = pr + other.pr;
    
    return result;
}

void Parameters::randInit() {
    static mt19937 gen(42); // Фиксируем seed
    // Определяем диапазоны для каждого параметра
    uniform_real_distribution<> dist_A1(0.0, 0.1);
    uniform_real_distribution<> dist_A0(0.0685, 0.1370);
    uniform_real_distribution<> dist_xi(0.7853, 1.57066667);
    uniform_real_distribution<> dist_pr(7.2853, 14.5706);
    uniform_real_distribution<> dist_q(2.0927, 4.1853);
    uniform_real_distribution<> dist_r0(1.9257, 3.8514);

    // Генерируем значения в соответствующих диапазонах
    A1 = dist_A1(gen);
    A0 = dist_A0(gen);
    xi = dist_xi(gen);
    pr = dist_pr(gen);
    q = dist_q(gen);
    r0 = dist_r0(gen);
}