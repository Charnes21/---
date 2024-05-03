#include <iostream>
#include <cmath>
#include <complex>
#include <locale.h>

using namespace std;

const int MAX_ROOTS = 100;
//Трансцендентное уравнения
double f11(double x)
{
    return(x * x - cos(5 * x));
}
double f111(double x)//Первая производная функции f11(x). 
{
    return(2 * x + 5 * sin(5 * x));
}
double f112(double x)  //Вторая производная функции f11(x). 
{
    return(2 + 25 * cos(5 * x));
}
double f11i(double x, double L)
{
    return(x + L * f11(x));
}

// Функция для вычисления значения многочлена с помощью схемы Горнера
double horner(double coefficients[], int n, double x) {
    double result = coefficients[0];


    for (int i = 1; i < n; ++i) {
        result = result * x + coefficients[i];
    }

    return result;
}

// Функция для деления многочлена на (x - a) и обновления его коэффициентов
void delPolynom(double coefficients[], int& n, double a) {
    if (horner(coefficients, n + 1, a) == 0) {
        // Применяем схему Горнера для деления многочлена на (x - a)
        double nextCoefficient = coefficients[0];
        coefficients[0] = nextCoefficient;
        for (int i = 1; i < n; ++i) {
            double currentCoefficient = coefficients[i];
            coefficients[i] = nextCoefficient;
            nextCoefficient = currentCoefficient + a * nextCoefficient; // вычисляем следующий коэффициент схемы Горнера
        }

        // Уменьшаем степень многочлена
        --n;
    }
}

// Алгебраические уравнение
double f1(double x) {
    return x * x * x * x - 10 * x * x * x + 35 * x * x - 50 * x + 24;
}

// Производная аглебр. уравн.
double proz(double x) {
    return 4 * pow(x, 3) - 30 * pow(x, 2) + 70 * x - 50;
}

// Функция, реализующая метод Ньютона
int newton(double st[], int numStPoints, double toch, double roots[]) {
    int num_roots = 0;
    for (int i = 0; i < numStPoints; ++i) {
        double x = st[i];
        int it = 0;
        bool found_root = false;
        while (fabs(f1(x)) > toch && num_roots < MAX_ROOTS) {
            x = x - f1(x) / proz(x);
            it++;
        }
        for (int j = 0; j < num_roots; ++j) {
            if (fabs(roots[j] - x) < 1e-10) {
                found_root = true;
                break;
            }
        }
        if (!found_root) {
            roots[num_roots++] = x;
            cout << "Количество итераций для начального приближения " << i + 1 << ": " << it << endl;
        }
    }
    return num_roots;
}

//вычисления разделенной разности первого порядка  
complex<double> razd(complex<double> x1, complex<double> x2, complex<double> f1, complex<double> f2) {
    return (f2 - f1) / (x2 - x1);
}

//вычисления разделенной разности второго порядка  
complex<double> razd2(complex<double> x0, complex<double> x1, complex<double> x2, complex<double> f0, complex<double> f1, complex<double> f2) {
    return (razd(x1, x2, f1, f2) - razd(x0, x1, f0, f1)) / (x0 - x2);
}

//вычисления значения полинома в точке  
complex<double> gorner(complex<double> P[], int N, complex<double> x0, int pr) {
    complex<double> y0 = 0;
    if (pr == 1) {
        for (int i = 0; i < N; ++i) {
            y0 = P[i] + x0 * y0;
        }
    }
    else if (pr == 2) {
        for (int i = 0; i < N; ++i) {
            y0 += pow(x0, N - i - 1) * P[i];
        }
    }
    return y0;
}

// функция понижения порядка полинома
void del_gor(complex<double> P[], int N, complex<double> alfa, complex<double> B[]) {
    B[0] = P[0];
    for (int i = 1; i < N; ++i) {
        B[i] = B[i - 1] * alfa + P[i];
    }
}

// методом парабол  
complex<double> parabol(complex<double> P[], int N, complex<double> x0, complex<double> x1, complex<double> x2) {
    complex<double> f0 = gorner(P, N, x0, 1);
    complex<double> f1 = gorner(P, N, x1, 1);
    complex<double> f2 = gorner(P, N, x2, 1);
    complex<double> x = x2;
    complex<double> eps = 1e-10;
    while (abs(gorner(P, N, x, 1)) > abs(eps)) {
        complex<double> u = razd(x0, x1, f0, f1) + razd(x0, x2, f0, f2) - razd(x1, x2, f1, f2);
        complex<double> t1 = u - sqrt(u * u - 4.0 * f0 * razd2(x0, x1, x2, f0, f1, f2));
        complex<double> t2 = u + sqrt(u * u - 4.0 * f0 * razd2(x0, x1, x2, f0, f1, f2));
        complex<double> t = (abs(t1) > abs(t2)) ? t1 : t2;
        x = x0 - 2.0 * gorner(P, N, x0, 1) / t;
        x2 = x1; x1 = x0; x0 = x;
        f0 = gorner(P, N, x0, 1);
        f1 = gorner(P, N, x1, 1);
        f2 = gorner(P, N, x2, 1);
    }
    return x;
}

//нахождения всех корней полинома  
void korni_polynom(complex<double> P[], int N, complex<double> x0, complex<double> x1, complex<double> x2) {
    int t = 0;
    cout << "Корни полинома:\n";
    while (N > 1) {
        complex<double> x = parabol(P, N, x0, x1, x2);
        cout << "x = " << x << endl;
        complex<double>* B = new complex<double>[N - 1];
        del_gor(P, N, x, B);
        for (int i = 0; i < N - 1; ++i) {
            P[i] = B[i];
        }
        --N;
        ++t;
    }
}



// Функция, заданная выражением (4.4)
double fi(double x, double L) {
    return x + L * (x * x - cos(5 * x));
}

// Функция, реализующая метод ложного положения
int Iteration(double* x, double L, double eps, double (*fi_)(double, double)) {
    int k = 0;
    double x0;
    do {
        x0 = *x;
        *x = fi_(*x, L);
        k++;
    } while (fabs(x0 - *x) >= eps);
    return k;
}

// Функция, реализующая метод половинного деления
int Dichotomy(double a, double b, double* c, double eps, double (*f_)(double)) {
    int k = 0;
    do {
        *c = (a + b) / 2;
        if (f_(*c) * f_(a) < 0) b = *c;
        else a = *c;
        k++;
    } while (fabs(a - b) >= eps);
    return k;
}

// Функция, реализующая метод хорд
int Chord(double a, double b, double* c, double eps, double (*f_)(double)) {
    int k = 0;
    do {
        *c = a - f_(a) / (f_(b) - f_(a)) * (b - a);
        if (f_(*c) * f_(a) > 0) a = *c;
        else b = *c;
        k++;
    } while (fabs(f_(*c)) >= eps);
    return k;
}

// Функция, реализующая метод касательных
int Tangent(double a, double b, double* c, double eps, double(*f_)(double), double(*f1_)(double), double(*f2_)(double)) {
    int k = 0;
    if (f_(a) * f2_(a) > 0) *c = a;
    else *c = b;
    do {
        *c = *c - f_(*c) / f1_(*c);
        k++;
    } while (fabs(f_(*c)) >= eps);
    return k;
}
// Функция ϕ(x), заданная по методу простой итерации
double phi(double x, double lambda) {
    return x + lambda * f11(x);
}

// Метод простой итерации для решения уравнения f(x) = 0
double simple_iteration(double A, double B, double lambda, double epsilon) {
    double x0 = (A + B) / 2; // Начальное приближение
    double x1 = phi(x0, lambda); // Первое приближение
    int iterations = 1; // Счетчик итераций

    // Проверяем условие сходимости на всем интервале [a, b]
    if (abs(1 + lambda * f111(A)) >= 1 || abs(1 + lambda * f111(B)) >= 1) {
        cout << "Метод не сходится на интервале [" << A << ", " << B << "]" << endl;
        return NAN;
    }

    // Итерационный процесс
    while (abs(x1 - x0) >= epsilon) {
        x0 = x1;
        x1 = phi(x0, lambda);
        iterations++;
    }

    cout << "Количество итераций: " << iterations << endl;
    return x1;
}

