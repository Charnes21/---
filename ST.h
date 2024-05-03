const int MAX_ROOTS = 100;
//Трансцендентное уравнения
double f11(double x);
double f111(double x);//Первая производная функции f11(x). 
double f112(double x); //Вторая производная функции f11(x). 

double f11i(double x, double L);

// Функция для вычисления значения многочлена с помощью схемы Горнера
double horner(double coefficients[], int n, double x);
// Функция для деления многочлена на (x - a) и обновления его коэффициентов
void delPolynom(double coefficients[], int& n, double a);

// Алгебраические уравнение
double f1(double x);

// Производная аглебр. уравн.
double proz(double x);

// Функция, реализующая метод Ньютона
int newton(double st[], int numStPoints, double toch, double roots[]);

//вычисления разделенной разности первого порядка  
complex<double> razd(complex<double> x1, complex<double> x2, complex<double> f1, complex<double> f2);

//вычисления разделенной разности второго порядка  
complex<double> razd2(complex<double> x0, complex<double> x1, complex<double> x2, complex<double> f0, complex<double> f1, complex<double> f2);

//вычисления значения полинома в точке  
complex<double> gorner(complex<double> P[], int N, complex<double> x0, int pr);

// функция понижения порядка полинома
void del_gor(complex<double> P[], int N, complex<double> alfa, complex<double> B[]);

// методом парабол  
complex<double> parabol(complex<double> P[], int N, complex<double> x0, complex<double> x1, complex<double> x2);

//нахождения всех корней полинома  
void korni_polynom(complex<double> P[], int N, complex<double> x0, complex<double> x1, complex<double> x2);


// Функция, заданная выражением (4.4)
double fi(double x, double L);
// Функция, реализующая метод ложного положения
int Iteration(double* x, double L, double eps, double (*fi_)(double, double));

// Функция, реализующая метод половинного деления
int Dichotomy(double a, double b, double* c, double eps, double (*f_)(double));

// Функция, реализующая метод хорд
int Chord(double a, double b, double* c, double eps, double (*f_)(double));

// Функция, реализующая метод касательных
int Tangent(double a, double b, double* c, double eps, double(*f_)(double), double(*f1_)(double), double(*f2_)(double));
// Функция ϕ(x), заданная по методу простой итерации
double phi(double x, double lambda);
// Метод простой итерации для решения уравнения f(x) = 0
double simple_iteration(double A, double B, double lambda, double epsilon);

