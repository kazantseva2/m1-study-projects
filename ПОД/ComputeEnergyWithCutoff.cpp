#include<iostream>
#include <cmath>  
#include <fstream>

// не считаем энергию дважды, "придуманные" атомы сами по себе

using namespace std;

const int N = 1; // сторона решетки
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
    // Для определения порядка атомов будем идти в следующем порядке:
    // сначала двигаемся по оси X до упора, 
    // затем делаем шаг по оси Y до тех пор пока не достигним границы,
    // затем сдвигаемся по оси Z.

    // Следовательно атомы в кубе будут обходиться в следующем порядке:
    // (0, 0, 0), (0.5a, 0.5a, 0), (0.5a, 0, 0.5a), (0, 0.5a, 0.5a).

    // После того, как все атомы упорядочены мы можем однозначно определить 
    // расстояние между i и j атомом.

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

Vector getNextAtom(Vector cur, double a) { // мб сделать метод, который не создает, а который меняет вектор
    if (cur.x == -1) {
        return Vector(0, 0, 0); //костыль
    }

    double border = N*a;
    if (cur.x + a < border) {
        return Vector(cur.x + a, cur.y, cur.z);
    }

    if (cur.x > border - a) {
        cur.x = 0;
    } else {
        cur.x = a/2;
    }

    if (cur.y + a/2 < border) {
        return Vector(cur.x, cur.y + a/2, cur.z);
    }

    cur.x = a/2;
    cur.y = 0;

    if (cur.z + a/2 < border) {
        return Vector(cur.x, cur.y, cur.z + a/2);
    }
    
    cout << "Error: fall outside" << endl;
    return Vector(0, 0, 0); //ужасно
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
    // cur.print();
    // cout << " <> ";
    // about.print();
    // cout << endl;

    if (about.x < 0 || about.x >= N * a || 
        about.y < 0 || about.y >= N * a || 
        about.z < 0 || about.z >= N * a)
        return false; // атом about не мог быть cur, тк он "придуманный"

    if (about.z != cur.z)
        return about.z < cur.z; // если координата z у текущего больше => был about для предыдущего, иначе не был

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

const double getDistanceBetweenAtoms(Vector &v1, Vector &v2) {
    if (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z) // сделать метод equals для Vector
        return 0;
        
    double rx = v1.x - v2.x; // сделать разность и покоординатное возведение в квадрат(или умножение пкоординатное)
    double ry = v1.y - v2.y;
    double rz = v1.z - v2.z;

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

double computeEnergyWithCutoff(double a) {
    double Eb, Er, E = 0;
    double cutoff = 1.7 * a, radCut = a; // если cutoff не меняется, то можно считать быстрее
    cutoff = a / sqrt(2); radCut = a / 2;
    // cout << "cutoff = " << cutoff << endl;
    
    double x_border = N * a, y_border = N * a, z_border = N * a; // зачем?
    double xi = 1.178;
    double r0 = 4.085 / sqrt(2);

    Vector v1 = Vector(-1, -1, -1);
    for (int i = 0; i < M * N * N * N; i++) {

        // if (i > 10)
        //     break;
        Eb = Er = 0;

        v1 = getNextAtom(v1, a);
        cout << "cur = ";
        v1.print();
        cout << endl;
        Vector v2 = Vector(v1);
        int n_around = 0;
        for (v2.x = v1.x - radCut; v2.x <= v1.x + radCut; v2.x += a/2) { 
            for (v2.y = v1.y - radCut; v2.y <= v1.y + radCut; v2.y += a/2) {
                for (v2.z = v1.z - radCut; v2.z <= v1.z + radCut; v2.z += a/2) {
                    // v1.print();
                    // cout << " ---> ";
                    // v2.print();

                    // считаем энергию 
                    if (!isAtom(v2, a) || alreadyCounted(v1, v2, a)) {
                        // cout << " continue" << endl;
                        continue;
                    }

                    double r = getDistanceBetweenAtoms(v1, v2); // TODO: сделать 'a' полем класса
                    if (r == 0) { // mb add cutoff
                        // cout << " r = " << r << endl;
                        continue;
                    }


                    
                    // cout << " r = " << r << endl;
                    n_around += 1;
                    // v1.print();
                    // cout << " ---> ";
                    // v2.print();
                    // cout << endl;

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
    
    return 0;
}