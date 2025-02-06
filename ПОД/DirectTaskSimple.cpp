#include<iostream>
#include <cmath>  
#include <fstream>

using namespace std;

const int N = 3; // сторона решетки
const int M = 4; // колличество базисных атомов 

struct Vector {
    double x, y, z;
    Vector(double x0, double y0, double z0) {
        x = x0;
        y = y0;
        z = z0;
    };

    void print() {
        cout << "(" << x << ", " << y  << ", "<< z <<")";
    }
};

Vector getCoordinates(int i, double a) {
    int nCube = i / 4; // номер куба в крисатлической решетке
    int nAtom = i % 4; // номер атома в кубе

    // по номеру куба определяем его координаты
    int zCube = nCube / (N*N);
    int yCube = (nCube % (N*N)) / N;
    int xCube = (nCube % (N*N)) % N;

    // определяем координаты атома
    // (0, 0, 0) --0, (0.5a, 0.5a, 0)--1, (0.5a, 0, 0.5a)--2, (0, 0.5a, 0.5a)--3.
    double x = xCube*a + 0.5*a*(nAtom == 1 or nAtom == 2);
    double y = yCube*a + 0.5*a*(nAtom == 1 or nAtom == 3);
    double z = zCube*a + 0.5*a*(nAtom == 2 or nAtom == 3);

    return Vector(x, y, z);
}

bool isInt(double a) {
    return a == floor(a);
}

const bool isAtom(Vector &v, double a) {
    double coefX = v.x / a;
    double coefY = v.y / a;
    double coefZ = v.z / a;

   // cout << "(" << coefX << ", " << coefY  << ", "<< coefZ <<")" << endl;
   // cout << "coef --> (" << isInt(coefX) << ", " << isInt(coefY)  << ", "<< isInt(coefZ) <<")" << endl;

    if (isInt(coefX) && isInt(coefY) && isInt(coefZ))
        return true;

    if (isInt(coefX) && !isInt(coefY) && !isInt(coefZ))
        return true;

    if (!isInt(coefX) && isInt(coefY) && !isInt(coefZ))
        return true;

    if (!isInt(coefX) && !isInt(coefY) && isInt(coefZ))
        return true;


    // cout << "FALSE" << endl;
    return false;
}

const bool alreadyCounted(Vector &cur, Vector &about, double a) { // mb error
    if (about.x < 0 || about.x >= N * a || 
        about.y < 0 || about.y >= N * a || 
        about.z < 0 || about.z >= N * a)
        return false;

    if (about.z < cur.z)
        return true;

    if (about.z > cur.z)
        return false;

    if (about.y < cur.y)
        return true;

    if (about.y > cur.y)
        return false;

    if (about.x < cur.x)
        return true;

    if (about.x > cur.x)
        return false;

    return true;
}


double cyclic_abs(double x, double a) {
    x = abs(x);
    if (x > N*a) {
        return N*a - x;
    }
    return x;
}

double getDistanceBetweenAtoms(int i, int j, double a) {
    if (i == j) 
        return 0;

    Vector v1 = getCoordinates(i, a); // NB: возможно вынести выше
    Vector v2 = getCoordinates(j, a);
        
    double rx = cyclic_abs(v1.x - v2.x, a);
    double ry = cyclic_abs(v1.y - v2.y, a);
    double rz = cyclic_abs(v1.z - v2.z, a);

    return sqrt(rx*rx + ry*ry + rz*rz);
}

const double getDistanceBetweenAtoms(Vector &v1, Vector &v2, double a) {
    if (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z) // сделать метод equals для Vector
        return 0;
        
    double rx = cyclic_abs(v1.x - v2.x, a);
    double ry = cyclic_abs(v1.y - v2.y, a);
    double rz = cyclic_abs(v1.z - v2.z, a);

    return sqrt(rx*rx + ry*ry + rz*rz);
}

double computeBindingEnergyForAtomIJ(double r, double r0) {
    double q = 3.139;
    double tmp = r/r0 - 1;
    return exp(-2 * q * tmp);
}

double computeRepulsionEnergyForAtomIJ(double r, double r0) {
    double A1 = 0.1028, A0 = 0.1028, p = 10.928;
    A0 = 0.1028, A1 = 0.000001, p = 10.928;
    double tmp = r/r0 - 1;
    return (A1 * (r - r0) + A0) * exp(-1 * p * tmp);
}

double computeEnergy(double a) {
    double Eb, Er, E = 0;
    double cutoff = 1.7 * a;
    double x_border = N * a, y_border = N * a, z_border = N * a;
    //double xi = 1.178;
    double xi = 1.178;
    double r0 = 4.085 / sqrt(2);
    //double r0 = a / sqrt(2);
    //double a = r0 * sqrt(2);



    for (int i = 0; i < M * N * N * N; i++) {
        Eb = Er = 0;
        for (int j = i; j < M * N * N * N; j++) {
            if (i == j) 
                continue;
            
            double r = getDistanceBetweenAtoms(i, j, a); // TODO: сделать 'a' полем класса
            cout << "dist(" << i << ", " << j << ") = " << r << endl;
            Eb += computeBindingEnergyForAtomIJ(r, r0);
            Er += computeRepulsionEnergyForAtomIJ(r, r0);
        } 
        E +=  Er - xi*sqrt(Eb);
    } 
    return E;
}

double computeEnergyWithCutoff(double a) {
    double Eb, Er, E = 0;
    double cutoff = 1.7 * a, radCut = a; // если cutoff не меняется, то можно считать быстрее
    cutoff = a / sqrt(2); radCut = a / 2;
    // cout << "cutoff = " << cutoff << endl;
    
    double x_border = N * a, y_border = N * a, z_border = N * a; // зачем?
    double xi = 1.178;
    double r0 = 4.085 / sqrt(2);



    for (int i = 0; i < M * N * N * N; i++) {
        if (i > 10)
            break;
        Eb = Er = 0;

        Vector v1 = getCoordinates(i, a);
        Vector v2 = Vector(v1);
        int n_around = 0;
        for (v2.x = v1.x - radCut; v2.x <= v1.x + radCut; v2.x += a/2) { 
            for (v2.y = v1.y - radCut; v2.y <= v1.y + radCut; v2.y += a/2) {
                for (v2.z = v1.z - radCut; v2.z <= v1.z + radCut; v2.z += a/2) {
                    // v1.print();
                    // cout << " ---> ";
                    // v2.print();

                    // x_j, y_j, z_j - какой-то постепенный перебор
                    // должна быть проверка: является ли точка с такой координатой атомом в решеьке( если это не решается обходом)

                    // считаем энергию 
                    if (!isAtom(v2, a) || alreadyCounted(v1, v2, a)) {
                        // cout << " continue" << endl;
                        continue;
                    }

                    double r = getDistanceBetweenAtoms(v1, v2, a); // TODO: сделать 'a' полем класса
                    if (r == 0) { // mb add cutoff
                        // cout << " r = " << r << endl;
                        continue;
                    }


                    
                    // cout << " r = " << r << endl;
                    n_around += 1;

                    Eb += computeBindingEnergyForAtomIJ(r, r0);
                    Er += computeRepulsionEnergyForAtomIJ(r, r0);
                }
            }
        }
        cout << n_around << endl;

        E +=  Er - xi*sqrt(Eb);
    } 
    return E;
}

int main()
{
    int P = N*N*N*M;
    double a = 1; //4.085;



    double y = computeEnergyWithCutoff(a) / P;
    cout << "a = " << a << " E = " << y << endl;

    // ofstream outputFile("data.txt");
    // if (outputFile.is_open()) {
    //     cout << "Увеличиваем" << endl;
    //     for (int i = 0; i < 100; i++) {
    //         a +=  0.01;
    //         double y = computeEnergyWithCutoff(a) / P;
    //         outputFile << a  << " " << y << "\n";
    //         cout << "a = " << a << " E = " << y << endl;
    //     }

    //     a = 4.085;
    //     cout << "Уменьшаем" << endl;
    //     for (int i = 0; i < 100; i++) {
    //         a -=  0.01;
    //         double y = computeEnergyWithCutoff(a) / P;
    //         outputFile << a  << " " << y << "\n";
    //         cout << "a = " << a << " E = " << y << endl;
    //     }

    //     a = 4.085;
    //     cout << "равномерно" << endl;
    //     for (int i = 0; i < 100; i++) {
    //         a +=  0.01 * i;
    //         double y = computeEnergyWithCutoff(a) / P;
    //         outputFile << a  << " " << y << "\n";
    //     }
    //     outputFile.close();
    // } else {
    //     cerr << "Не удалось открыть файл для записи!" << endl;
    // }
    
    return 0;
}