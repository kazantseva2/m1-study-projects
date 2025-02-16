#include "definitions.h"


Vector::Vector(double x0, double y0, double z0) : x(x0), y(y0), z(z0) {}

Vector::Vector() : x(0), y(0), z(0) {}

Vector Vector::add(double a, double b, double c) const {
    return Vector(x + a, y + b, z + c);
}

bool Vector::operator==(Vector other) const {
    return this->x == other.x && this->y == other.y && this->z == other.z;
}

Vector Vector::operator-(Vector other) const {
    return Vector(this->x - other.x, this->y - other.y, this->z - other.z);
}

Vector Vector::operator*(Vector other) const {
    return Vector(this->x * other.x, this->y * other.y, this->z * other.z);
}

double Vector::length() const {
    return sqrt(x * x + y * y + z * z);
}

double Vector::mod_x(double x, double n) const {
    x = abs(x);
    return min(x, n - x);
}

Vector Vector::mod(double n, bool vacuum) const {
    double mx = mod_x(this->x, n);
    double my = mod_x(this->y, n);
    double mz;
    if (vacuum) {
        mz = abs(this->z);
    } else {
        mz = mod_x(this->z, n);
    }

    return Vector(mx, my, mz);
}

Vector Vector::operator*(const vector<vector<double>>& matrix) const {
    return Vector(
        x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2],
        x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2],
        x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2]
    );
}