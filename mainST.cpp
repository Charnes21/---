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
    cout << "�������� ����� ������� ���������:" << endl;
    cout << "1. ����� ��������" << endl;
    cout << "2. ����� �������" << endl;
    cout << "3. ������ ����������� �������, ����, �����������, ������� ��������� � ����� ������� �������� " << endl;
    cout << "��� �����: ";
    cin >> choice;

    switch (choice) {
    case 1: {
        int size;
        cout << "������� ������� ��������: ";
        cin >> size;
        ++size;

        complex<double>* H = new complex<double>[size];
        cout << "������� ������������ �������� (������� �� �������� �����):" << endl;
        for (int i = 0; i < size; ++i) {
            cout << "����������� ��� x^" << (size - i - 1) << ": ";
            cin >> H[i];
        }

        complex<double> x0, x1, x2;
        cout << "������� ��������� ������������� x0: " << endl;
        cout << "���� ����� ������ � ������� ����������� ����� (2, 1), ��� 1 - ������ �����, ���� x����� ������ ������������ �����, ������ ����� ������ = 0 " << endl;

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

        cout << "������� ����������� �������� x: ";
        cin >> minX;
        cout << "������� ������������ �������� x: ";
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
        cout << "������� ���������� ������, ������� �� �������� �����: ";
        cin >> pribliz;

        double stPoints[MAX_ROOTS];
        cout << "������� " << pribliz << " ��������� �����������: ";
        for (int i = 0; i < pribliz; ++i) {
            cin >> stPoints[i];
        }

        double toch;
        cout << "������� ��������: ";
        cin >> toch;

        double kor[MAX_ROOTS];
        int num_roots = newton(stPoints, pribliz, toch, kor);

        cout << "���������� ��������� ������ ���������: " << num_roots << endl;
        cout << "����� ���������: ";
        for (int i = 0; i < num_roots; ++i) {
            cout << kor[i] << " ";
        }
        cout << endl;

        break;
    }
    case 3: {
        double minX, maxX;

        cout << "������� ����������� �������� x: ";
        cin >> minX;
        cout << "������� ������������ �������� x: ";
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
        double ep = 0.001;  //�������� ����������. 
        int K;
        cout << "��������� x * x - cos(5 * x)" << endl;
        cout << "������� ������ �������\n";
        cout << "a="; cin >> A;  //�������� �������� �����. 
        cout << "������� ����� �������\n";
        cout << "b="; cin >> B;
        cout << "������� ��������� x^2-cos(5*x)=0." << endl;
        cout << "����� ����������� �������:" << endl;
        K = Dichotomy(A, B, &X, ep, f11);
        cout << "�������  ���������  x=" << X;
        cout << ", ���������� �������� k=" << K << endl;
        cout << "�����  ����:" << endl;
        K = Chord(A, B, &X, ep, f11);
        cout << "�������  ���������  x=" << X;
        cout << ", ���������� �������� k=" << K << endl;
        cout << " ����� ������� ���������:" << endl;
        X = A;
        cout << "������� ��������� ������������� �����" << endl;
        cout << "L="; cin >> P;
        K = Iteration(&X, P, ep, f11i);
        cout << "�������  ���������  x=" << X;
        cout << ", ���������� �������� k=" << K << endl;
        double lambda;
        double epsilon;

        cout << "������� �������� ��������� lambda: ";
        cin >> lambda;
        cout << "������� �������� �������� epsilon: ";
        cin >> epsilon;

        double root = simple_iteration(A, B, lambda, epsilon);
        if (!isnan(root)) {
            cout << "������ ���������: " << root << endl;
        }
        else {
            cout << "�� ������� ����� ������ �� ��������� ���������." << endl;
        }
        break;
    }
    default:
        cout << "������: �������� ����� ������." << endl;
        break;
    }

    return 0;
}

