#include "classes.h"
#include <iostream>
#include <cmath>


using namespace std;

Vector::Vector(double x0, double y0, double z0) : x(x0), y(y0), z(z0) {}

Vector::Vector() : x(0), y(0), z(0) {}

void Vector::print() const {
    cout << "(" << x << ", " << y << ", " << z << ")";
}

void Vector::modf() {
    x = abs(x - static_cast<int>(x));
    y = abs(y - static_cast<int>(y));
    z = abs(z - static_cast<int>(z));
}

bool Vector::equals(const Vector& v) const {
    return x == v.x && y == v.y && z == v.z;
}

double Vector::length() const {
    return sqrt(x * x + y * y + z * z);
}

Vector Vector::add(double a, double b, double c) const {
    return Vector(x + a, y + b, z + c);
}

Vector Vector::multiplyByVector(const Vector& coefs) const {
    return Vector(x * coefs.x, y * coefs.y, z * coefs.z);
}

Vector Vector::multiplyByMatrix(const vector<vector<double>>& matrix) const {
    return Vector(
        x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2],
        x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2],
        x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2]
    );
}

bool Vector::operator==(const Vector& other) const {
    return x == other.x && y == other.y && z == other.z;
}

Vector Vector::operator-(const Vector& other) const {
    return Vector(x - other.x,  y - other.y, z - other.z);
}
