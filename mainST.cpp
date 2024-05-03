#include <iostream>
#include <cmath>
#include <complex>
#include <locale.h>
#include "GnuP.h"
#include "ST.h"

using namespace std;

int main() {
    setlocale(LC_ALL, "Russian");

    int choice;
    cout << "Выберите метод решения уравнения:" << endl;
    cout << "1. Метод Параболы" << endl;
    cout << "2. Метод Ньютона" << endl;
    cout << "3. Методы половинного деления, хорд, касательных, ложного положения и метод простой итерации " << endl;
    cout << "Ваш выбор: ";
    cin >> choice;

    switch (choice) {
    case 1: {
        int size;
        cout << "Введите степень полинома: ";
        cin >> size;
        ++size;

        complex<double>* H = new complex<double>[size];
        cout << "Введите коэффициенты полинома (начиная со старшего члена):" << endl;
        for (int i = 0; i < size; ++i) {
            cout << "Коэффициент для x^" << (size - i - 1) << ": ";
            cin >> H[i];
        }

        complex<double> x0, x1, x2;
        cout << "Введите начальное предположение x0: " << endl;
        cout << "Ввод нужно ввести в формате комлексного числа (2, 1), где 1 - мнимая часть, если xотите ввести вещественное число, мнимую часть запиши = 0 " << endl;

        double real, imag;
        cin >> real >> imag;
        x0 = complex<double>(real, imag);
        x1 = x0 * 0.85;
        x2 = x0 * 1.15;

        korni_polynom(H, size, x0, x1, x2);
        delete[] H;

        break;
    }
    case 2: {
        double minX, maxX;

        cout << "Введите минимальное значение x: ";
        cin >> minX;
        cout << "Введите максимальное значение x: ";
        cin >> maxX;

        int numPoints = static_cast<int>((maxX - minX) / 0.1) + 1;

        double* x = new double[numPoints];
        double* y = new double[numPoints];

        for (int i = 0; i < numPoints; ++i)
        {
            x[i] = minX + i * 0.1;
            y[i] = f1(x[i]);
        }

        GnuP p;
        p.plotArray(numPoints, x, y);
        p.plot();

        delete[] x;
        delete[] y;
        //double st;
        int pribliz;
        cout << "Введите количество корней, которые вы ожидаете найти: ";
        cin >> pribliz;

        double stPoints[MAX_ROOTS];
        cout << "Введите " << pribliz << " начальных приближений: ";
        for (int i = 0; i < pribliz; ++i) {
            cin >> stPoints[i];
        }

        double toch;
        cout << "Введите точность: ";
        cin >> toch;

        double kor[MAX_ROOTS];
        int num_roots = newton(stPoints, pribliz, toch, kor);

        cout << "Количество найденных корней уравнения: " << num_roots << endl;
        cout << "Корни уравнения: ";
        for (int i = 0; i < num_roots; ++i) {
            cout << kor[i] << " ";
        }
        cout << endl;

        break;
    }
    case 3: {
        double minX, maxX;

        cout << "Введите минимальное значение x: ";
        cin >> minX;
        cout << "Введите максимальное значение x: ";
        cin >> maxX;

        int numPoints = static_cast<int>((maxX - minX) / 0.1) + 1;

        double* x = new double[numPoints];
        double* y = new double[numPoints];

        for (int i = 0; i < numPoints; ++i)
        {
            x[i] = minX + i * 0.1;
            y[i] = f11(x[i]);
        }

        GnuP p;
        p.plotArray(numPoints, x, y);
        p.plot();

        delete[] x;
        delete[] y;
        double A, B, X, P;
        double ep = 0.001;  //Точность вычислений. 
        int K;
        cout << "уравнение x * x - cos(5 * x)" << endl;
        cout << "Введите начало отрезка\n";
        cout << "a="; cin >> A;  //Интервал изоляции корня. 
        cout << "Введите конец отрезка\n";
        cout << "b="; cin >> B;
        cout << "Решение уравнения x^2-cos(5*x)=0." << endl;
        cout << "Метод половинного деления:" << endl;
        K = Dichotomy(A, B, &X, ep, f11);
        cout << "Решение  уравнения  x=" << X;
        cout << ", количество итераций k=" << K << endl;
        cout << "Метод  хорд:" << endl;
        K = Chord(A, B, &X, ep, f11);
        cout << "Решение  уравнения  x=" << X;
        cout << ", количество итераций k=" << K << endl;
        cout << " Метод ложного положения:" << endl;
        X = A;
        cout << "Введите начальное предположение корня" << endl;
        cout << "L="; cin >> P;
        K = Iteration(&X, P, ep, f11i);
        cout << "Решение  уравнения  x=" << X;
        cout << ", количество итераций k=" << K << endl;
        double lambda;
        double epsilon;

        cout << "Введите значение параметра lambda: ";
        cin >> lambda;
        cout << "Введите значение точности epsilon: ";
        cin >> epsilon;

        double root = simple_iteration(A, B, lambda, epsilon);
        if (!isnan(root)) {
            cout << "Корень уравнения: " << root << endl;
        }
        else {
            cout << "Не удалось найти корень на указанном интервале." << endl;
        }
        break;
    }
    default:
        cout << "Ошибка: неверный выбор метода." << endl;
        break;
    }

    return 0;
}

