#include<iostream>
#include <cmath>  
#include <fstream>

using namespace std;

const int N = 3; // сторона решетки
const int M = 4; // колличество базисных атомов 
const double E_coh = -2.96; // TODO: придумать, где задавать (тут для серебра)

struct Vector {
    double x, y, z;
    Vector(double x0, double y0, double z0) {
        x = x0;
        y = y0;
        z = z0;
    };

    Vector& operator-= (const Vector& o);
    Vector& operator*= (double o);
};

Vector operator- (const Vector& lhs, const Vector& rhs);
Vector operator* (const Vector& lhs, double rhs);
Vector operator* (double lhs, const Vector& rhs);
 
// struct Cube {
//     Vector atom000, atom110, atom011, atom101;  // basic atoms: (0, 0, 0), (0.5a, 0.5a, 0), 
//                                                             // (0, 0.5a, 0.5a), (0.5a, 0, 0.5a)

// };

struct GridModel {
   // Cube grid[N][N][N];
    double A0, A1, r0, q, p, xi; // параметры модели
    double a; // параметр решетки
    double E_coh; // когезионная энергия 

    // Для определения порядка атомов будем идти в следующем порядке:
    // сначала двигаемся по оси X до упора, 
    // затем делаем шаг по оси Y до тех пор пока не достигним границы,
    // затем сдвигаемся по оси Z.

    // Следовательно атомы в кубе будут обходиться в следующем порядке:
    // (0, 0, 0), (0.5a, 0.5a, 0), (0.5a, 0, 0.5a), (0, 0.5a, 0.5a).

    // После того, как все атомы упорядочены мы можем однозначно определить 
    // расстояние между i и j атомом.



    GridModel();

    const double computeBindingEnergyForAtomIJ(double r);
    const double computeRepulsionEnergyForAtomIJ(double r);
    double computeEnergy(); // рассчет полной энергии кристаллической ячейки

    const double cyclic_abs(double x); // NB: подумать над тем, чтоб вынести параметр а
    const Vector getCoordinates(int i);
    const double getDistanceBetweenAtoms(int i, int j);

};

GridModel::GridModel() {
    A0 = A1 = r0 = q = p = xi = 1;
    a = 1;
}


const Vector GridModel::getCoordinates(int i) {
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

const double GridModel::cyclic_abs(double x) {
    x = abs(x);
    if (x > N*a) {
        return N*a - x;
    }
    return x;
}
    
const double GridModel::getDistanceBetweenAtoms(int i, int j) {
    if (i == j) 
        return 0;

    Vector v1 = getCoordinates(i); // NB: возможно вынести выше
    Vector v2 = getCoordinates(j);
        
    double rx = cyclic_abs(v1.x - v2.x);
    double ry = cyclic_abs(v1.y - v2.y);
    double rz = cyclic_abs(v1.z - v2.z);

    return sqrt(rx*rx + ry*ry + rz*rz);
}

const double GridModel::computeBindingEnergyForAtomIJ(double r) {
    double tmp = r/r0 - 1;
    return exp(-2 * q * tmp);
}

double GridModel::computeRepulsionEnergyForAtomIJ(double r) {
    double tmp = r/r0 - 1;
    return (A1 * (r - r0) + A0) * exp(-1 * p * tmp);
}


double GridModel::computeEnergy() {
   // cout << "Calculate Energy" << endl;
    double Eb, Er, E = 0;
    Vector v1(0, 0, 0), v2(0, 0, 0);

    for (int i = 0; i < M * N * N * N; i++) {
        Eb = Er = 0;
        for (int j = 0; j < M * N * N * N; j++) {
            if (i == j) 
                continue;
            
            double r = getDistanceBetweenAtoms(i, j); 
            Eb += computeBindingEnergyForAtomIJ(r);
            Er += computeRepulsionEnergyForAtomIJ(r);
        } 
        E +=  Er - xi*sqrt(Eb);
    } 
    //cout << "Eb = " << Eb << ", Er = " << Er << endl;
    return E;
}


int main()
{
    int P = N*N*N*M;
    double r0 = 0.4;
    GridModel grid;

    // создаем решетку с каими-то параметрами, начнем с установки r0
    // считаем энергию 
    // пусть будет функция E(r0)
    // grid.r0 = r0;
    // cout << "a = " << grid.r0 << ", E = " << grid.computeEnergy() / P << endl;
    // for (int i = 0; i < 50; i++) {
    //     grid.r0 += grid.r0 * 0.01;
    //     cout << "a = " << grid.r0 << ", E = " << grid.computeEnergy() / P << endl;
    // }

    // cout << "уменьшаем" << endl;
    // grid.a = r0;
    // for (int i = 0; i < 30; i++) {
    //     grid.a -= grid.a * 0.01;
    //     cout << "a = " << grid.a << ", E = " << grid.computeEnergy() / P << endl;
    // }


    // cout << endl;
    // cout << "равномерно" << endl;
    // grid.a = r0;
    // for (int i = 0; i < 100; i++) {
    //     grid.a +=  0.01 * i;
    //     cout << "a = " << grid.a << ", E = " << grid.computeEnergy() / P << endl;
    // }

    // cout << "уменьшаем" << endl;
    // grid.a = r0;
    // for (int i = 0; i < 100; i++) {
    //     grid.a -= grid.a * 0.01 * i;
    //     cout << "a = " << grid.a << ", E = " << grid.computeEnergy() / P << endl;
    // }

    ofstream outputFile("data.txt");
    if (outputFile.is_open()) {
        cout << "равномерно" << endl;
        grid.r0 = r0;
        for (int i = 0; i < 100; i++) {
            grid.r0 +=  0.01 * i;
            double y = grid.computeEnergy() / P;
            outputFile << grid.r0  << " " << y << "\n";
        }
         cout << "уменьшаем" << endl;
        grid.a = r0;
        for (int i = 0; i < 100; i++) {
            grid.r0 -= grid.r * 0.01 * i;
            double y = grid.computeEnergy() / P;
            outputFile << grid.a  << " " << y << "\n";
        }
        outputFile.close();
    } else {
        std::cerr << "Не удалось открыть файл для записи!" << std::endl;
    }
    
    return 0;
}