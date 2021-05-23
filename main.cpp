#include <iostream>
#include <math.h>


using namespace std;

#define N 100
#define ALPHA 1E-8
#define BETA 1

void mygen(double **pDouble, double **pDouble1, int n, double alpha, double beta, int i, int i1, int i2, int i3);
void Q_matrix(double **pDouble, int n, int schema);
double matr_inf_norm(double **pDouble, int n);
void matr_mul(double **pDouble, double **pDouble1, double **pDouble2, int n);
double a(double x) {
    if (x >= 0)
        return x;
    else {
        x = x * (-1);
        return x;
    }
}

int maxStolb(int size, double **arr, int i) //поиск номера строки содержащей максимальный элемент
{
    int num_str = i;
    double max = a(arr[i][i]);
    for (int s = i; s < size; s++)
        if (a(arr[s][i]) > max) {
            num_str = s;
            max = a(arr[s][i]);
        }
    if (max != 0)
        return num_str;
    else { return -1; }
}

void perestStrok(int size, double **arr, int i, int j)//перестановка i-ой и j-ой строк
{
    double d;
    for (int s = 0; s < size; s++) {
        d = arr[i][s];
        arr[i][s] = arr[j][s];
        arr[j][s] = d;
    }
    return;
}

double *Gordan(int size, double **arr, double *f,double *xp) {
    int stl;
    for (int i = 0; i < size; i++) {
        stl = maxStolb(size, arr, i);
        if (stl < 0)
            return NULL;
        if (i != stl)
            perestStrok(size, arr, i, stl);
        double d = f[i];
        f[i] = f[stl];
        f[stl] = d;
        for (int k = 0; k < size; k++) {
            if (k != i) {
                double a = arr[k][i];
                for (int m = 0; m < size; m++)
                    arr[k][m] = arr[k][m] - arr[i][m] / arr[i][i] * a;
                f[k] = f[k] - f[i] / arr[i][i] * a;
            }
        }
    }
    for (int i = 0; i < size; i++)
        xp[i] = f[i] / arr[i][i];
    return xp;
}

double cube_norm(double *z, int n)
{
    double max = fabs(z[0]); //модуль числа
    for (int i = 1; i < n; i++){
        if (max < fabs(z[i])) max = fabs(z[i]);
    }
    return max;
}

//считает a*x(тр)-f
double *neviazka(double **a, double *xp, double *f, int n){
    double *res = new double [n];
    for (int i = 0; i < n; i++){
        res[i] = -f[i];
        for (int j = 0; j < n; j++) res[i] = res[i] + a[i][j] * xp[j];
    }
    return res;
}

//умножение матрицы на вектор
double *mul(double **&a, double *x, int n){
    double *res = new double[n];
    for (int i = 0; i < n; i++){
        res[i] = 0;
        for (int j = 0; j < n; j++) res[i] = res[i] + a[i][j] * x[j];
    }
    return res;
}

int main()
{
    int n = N;
    cout << " n = " << n << endl;

    //определить матрицу
    double **a = new double *[n];
    for (int i = 0; i < n; i++) a[i] = new double[n];
    double **a_inv = new double*[n];
    for (int i = 0; i < n; i++) a_inv[i] = new double[n];
    double alpha = ALPHA;
    double beta = BETA;

  // 	mygen ( a, a_inv, n, alpha, beta, 1, 2, 0, 1 ); // симметричная
    //mygen ( a, a_inv, n, alpha, beta, 1, 2, 1, 1 ); //простой структуры
    mygen(a, a_inv, n, alpha, beta, 0, 0, 2, 1); //жорданова клетка
    cout << endl;
    delete[] a_inv;

    //назначаем точное решение
    double *xt = new double[n];//точное решение
    xt[0] = -1.;
    xt[1] = 0.;
    xt[2] = sin(2.);
    for (int i = 3; i < n; i++)
        xt[i] = 10.;



    //правая часть
    double *f = mul(a, xt, n);

//выделяем память под приблеженное решение
    double *xp = new double[n];//приблеженное решение


    //определить матрицу 2
    double **a2 = new double *[n];
    for (int i = 0; i < n; i++) a2[i] = new double[n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a2[i][j] = a[i][j];


    //правая часть 2
    double *f2 = new double[n];
    for (int i = 0; i < n; i++){
        f2[i] = f[i];
    }

    //вызов решателя
    Gordan(n,a2,f2,xp);
    for (int i = 0; i < n; i++)
        delete[] a2[i];
    delete []a2;
    delete []f2;




        double *z = new double[n];

        for (int i = 0; i < n; i++) {
            z[i] = xp[i] - xt[i];
        }

        double z_norm = cube_norm(z, n);
        cout << "||z|| = " << z_norm << endl;

        //dzeta = ||z||/||xt|| ->на экран
        //r = A*xp (x с волной) - f
        //||r|| -> на экран
        //p = ||r||/||f|| - относительная норма

        double s = z_norm / cube_norm(xt, n);
        cout << "dzeta = " << s << endl;


        double *r = neviazka(a, xp, f , n);

        cout << "||r|| = " << cube_norm(r, n);
        double p = cube_norm(r, n) / cube_norm(f, n);
        cout << endl;

        cout << " p = " << p << endl;
        cout << endl;
        delete[]z;
        delete[]r;


    for (int i = 0; i < n; i++)
        delete[] a[i];
    delete[]a;

    delete[]f;
    delete[]xt;
    delete[]xp;

    return 0;
}

void mygen(double **a, double **a_inv, int n, double alpha, double beta, int sign_law, int lambda_law, int variant, int schema)
{
    int i, j, k;


    cout << "   M A T R I X  G E N.  " << endl;

    cout << "              N = " << n << endl;
    cout << " | lambda_min | = " << alpha << endl;
    cout << " | lambda_max | = " << beta << endl;

    double *lambda = new double[n];

    // распределение знаков
    cout << " sign_law = " << sign_law << endl;

    double *sign = new double[n];
    for (i = 0; i<n; i++) sign[i] = 1.;

    switch (sign_law)
    {
        case -1:
            for (i = 0; i<n; i++) sign[i] = -1.;
            break;

        case 0:
            sign[0] = 1.;
            for (i = 1; i<n; i++) sign[i] = -sign[i - 1];
            break;

            //другие законы распределения знаков
            // ...

    }

    //распределение собственнных чисел
    cout << " lambda_law = " << lambda_law << endl;

    double *kappa = new double[n];
    for (i = 0; i<n; i++) kappa[i] = (double)i / double(n - 1);
    switch (lambda_law)
    {
        case 1:
            cout << " kappa = sqrt( ) " << endl;
            for (i = 0; i<n; i++) kappa[i] = sqrt(kappa[i]);
            break;

        case 2:
            cout << " kappa = sin( ) " << endl;
            double pi_half = acos(-1.)*0.5;
            for (i = 0; i<n; i++) kappa[i] = sin(pi_half*kappa[i]);
            break;

    }


    double *J = new double[n];
    for (i = 0; i<n; i++) J[i] = sign[i] * ((1. - kappa[i])*alpha + kappa[i] * beta);

    double *J_inv = new double[n];
    for (i = 0; i<n; i++) J_inv[i] = 1. / J[i];

    double **Q = new double *[n];
    for (i = 0; i<n; i++) Q[i] = new double[n];

    double aa[3];


    cout << " variant = " << variant << endl;

    switch (variant)
    {
        case 0: //симметричная матрица
            cout << " simmetric matrix:" << endl;
            cout << " schema = " << schema << endl;
            switch (schema)
            {
                case 1:
                    Q_matrix(Q, n, schema);

                    for (a[0][0] = 0., k = 0; k<n; k++) a[0][0] += Q[0][k] * J[k] * Q[0][k];
                    for (j = 1; j<n; j++)
                    {
                        for (a[0][j] = 0., k = j - 1; k<n; k++) a[0][j] += Q[0][k] * J[k] * Q[j][k];
                        a[j][0] = a[0][j];
                    }
                    for (i = 1; i<n; i++)
                    {
                        for (a[i][i] = 0., k = i - 1; k<n; k++) a[i][i] += Q[i][k] * J[k] * Q[i][k];
                        for (j = i + 1; j<n; j++)
                        {
                            for (a[i][j] = 0., k = j - 1; k<n; k++) a[i][j] += Q[i][k] * J[k] * Q[j][k];
                            a[j][i] = a[i][j];
                        }
                    }

                    //_______
                    for (a_inv[0][0] = 0., k = 0; k<n; k++) a_inv[0][0] += Q[0][k] * J_inv[k] * Q[0][k];
                    for (j = 1; j<n; j++)
                    {
                        for (a_inv[0][j] = 0., k = j - 1; k<n; k++) a_inv[0][j] += Q[0][k] * J_inv[k] * Q[j][k];
                        a_inv[j][0] = a_inv[0][j];
                    }
                    for (i = 1; i<n; i++)
                    {
                        for (a_inv[i][i] = 0., k = i - 1; k<n; k++) a_inv[i][i] += Q[i][k] * J_inv[k] * Q[i][k];
                        for (j = i + 1; j<n; j++)
                        {
                            for (a_inv[i][j] = 0., k = j - 1; k<n; k++) a_inv[i][j] += Q[i][k] * J_inv[k] * Q[j][k];
                            a_inv[j][i] = a_inv[i][j];
                        }
                    }
                    break;

            }//schema
            break;

        case 1: //матрица простой структуры
            cout << " simple structure matrix:" << endl;
            cout << " schema = " << schema << endl;
            switch (schema)
            {
                case 1:
                    //TJ
                    //первая строка
                    a[0][0] = J[0];
                    a[0][1] = -J[1];
                    //			for(j=2; j<n; j++ ) a[0][j] = 0.;
                    //до последней
                    for (i = 1; i<n - 1; i++)
                    {
                        //				for(j=0; j<i-1; j++ ) a[i][j] = 0.;
                        a[i][i - 1] = -J[i - 1];
                        a[i][i] = J[i] + J[i];
                        a[i][i + 1] = -J[i + 1];
                        //				for(j=i+2; j<n; j++ ) a[i][j] = 0.;
                    }
                    //последняя (n-1)
                    //			for(j=0; j<n-2; j++ ) a[n-1][j] = 0.;
                    a[n - 1][n - 2] = -J[n - 2];
                    a[n - 1][n - 1] = J[n - 1] + J[n - 1];

                    //(TJ)T^{-1}
                    //первая строка
                    aa[1] = a[0][0];  aa[2] = a[0][1];
                    a[0][0] = aa[1] * (double)n + aa[2] * (double)(n - 1);
                    double s = aa[1] + aa[2];
                    for (j = 1; j<n; j++) a[0][j] = s*(double)(n - j);
                    //до последней
                    for (i = 1; i<n - 1; i++)
                    {
                        aa[0] = a[i][i - 1];  aa[1] = a[i][i];  aa[2] = a[i][i + 1];
                        for (j = 0; j<i; j++) a[i][j] = aa[0] * (double)(n - i + 1) + aa[1] * (double)(n - i) + aa[2] * (double)(n - i - 1);
                        s = aa[0] + aa[1];
                        a[i][i] = s*(double)(n - i) + aa[2] * (double)(n - i - 1);
                        s += aa[2];
                        for (j = i + 1; j<n; j++) a[i][j] = s*(double)(n - j);
                    }
                    //последняя (n-1)
                    aa[0] = a[n - 1][n - 2];  aa[1] = a[n - 1][n - 1];
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j<n - 1; j++) a[n - 1][j] = s;
                    a[n - 1][n - 1] = aa[0] + aa[1];
                    //_______

                    //TJ^{-1}
                    //первая строка
                    a_inv[0][0] = J_inv[0];
                    a_inv[0][1] = -J_inv[1];
                    //до последней
                    for (i = 1; i<n - 1; i++)
                    {
                        a_inv[i][i - 1] = -J_inv[i - 1];
                        a_inv[i][i] = J_inv[i] + J_inv[i];
                        a_inv[i][i + 1] = -J_inv[i + 1];
                    }
                    //последняя (n-1)
                    a_inv[n - 1][n - 2] = -J_inv[n - 2];
                    a_inv[n - 1][n - 1] = J_inv[n - 1] + J_inv[n - 1];

                    //(TJ^{-1})T^{-1}
                    //первая строка
                    aa[1] = a_inv[0][0];  aa[2] = a_inv[0][1];
                    a_inv[0][0] = aa[1] * (double)n + aa[2] * (double)(n - 1);
                    s = aa[1] + aa[2];
                    for (j = 1; j<n; j++) a_inv[0][j] = s*(double)(n - j);
                    //до последней
                    for (i = 1; i<n - 1; i++)
                    {
                        aa[0] = a_inv[i][i - 1];  aa[1] = a_inv[i][i];  aa[2] = a_inv[i][i + 1];
                        for (j = 0; j<i; j++) a_inv[i][j] = aa[0] * (double)(n - i + 1) + aa[1] * (double)(n - i) + aa[2] * (double)(n - i - 1);
                        s = aa[0] + aa[1];
                        a_inv[i][i] = s*(double)(n - i) + aa[2] * (double)(n - i - 1);
                        s += aa[2];
                        for (j = i + 1; j<n; j++) a_inv[i][j] = s*(double)(n - j);
                    }
                    //последняя (n-1)
                    aa[0] = a_inv[n - 1][n - 2];  aa[1] = a_inv[n - 1][n - 1];
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j<n - 1; j++) a_inv[n - 1][j] = s;
                    a_inv[n - 1][n - 1] = aa[0] + aa[1];
                    break;

            }//schema
            break;

        case 2: //одна жорданова клетка 2x2 при минимальном с.з.
            cout << " J_2 type matrix:" << endl;
            cout << " schema = " << schema << endl;

            switch (schema)
            {
                case 1:
                    //TJ
                    //первая строка
                    a[0][0] = J[0];
                    a[0][1] = 1. - J[0];
                    //вторая строка
                    a[1][0] = -J[0];
                    a[1][1] = -1. + J[0] + J[0];
                    a[1][2] = -J[2];
                    //третья строка
                    a[2][1] = -J[0];
                    a[2][2] = J[2] + J[2];
                    if (n>3) a[2][3] = -J[3];
                    //до последней
                    for (i = 3; i<n - 1; i++)
                    {
                        a[i][i - 1] = -J[i - 1];
                        a[i][i] = J[i] + J[i];
                        a[i][i + 1] = -J[i + 1];
                    }
                    //последняя (n-1)
                    if (n>3)
                    {
                        a[n - 1][n - 2] = -J[n - 2];
                        a[n - 1][n - 1] = J[n - 1] + J[n - 1];
                    }

                    //(TJ)T^{-1}
                    //первая строка
                    aa[1] = a[0][0];  aa[2] = a[0][1];
                    a[0][0] = aa[1] * (double)n + aa[2] * (double)(n - 1);
                    double s = aa[1] + aa[2];
                    for (j = 1; j<n; j++) a[0][j] = s*(double)(n - j);
                    //до последней
                    for (i = 1; i<n - 1; i++)
                    {
                        aa[0] = a[i][i - 1];  aa[1] = a[i][i];  aa[2] = a[i][i + 1];
                        for (j = 0; j<i; j++) a[i][j] = aa[0] * (double)(n - i + 1) + aa[1] * (double)(n - i) + aa[2] * (double)(n - i - 1);
                        s = aa[0] + aa[1];
                        a[i][i] = s*(double)(n - i) + aa[2] * (double)(n - i - 1);
                        s += aa[2];
                        for (j = i + 1; j<n; j++) a[i][j] = s*(double)(n - j);
                    }
                    //последняя (n-1)
                    aa[0] = a[n - 1][n - 2];  aa[1] = a[n - 1][n - 1];
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j<n - 1; j++) a[n - 1][j] = s;
                    a[n - 1][n - 1] = aa[0] + aa[1];
                    //_______
                    //TJ^{-1}
                    //первая строка
                    a_inv[0][0] = J_inv[0];
                    a_inv[0][1] = -J_inv[0] * J_inv[0] - J_inv[0];
                    //вторая строка
                    a_inv[1][0] = -J_inv[0];
                    a_inv[1][1] = J_inv[0] * J_inv[0] + J_inv[0] + J_inv[0];
                    a_inv[1][2] = -J_inv[2];
                    //третья строка
                    a_inv[2][1] = -J_inv[0];
                    a_inv[2][2] = J_inv[2] + J_inv[2];
                    if (n>3) a_inv[2][3] = -J_inv[3];
                    //до последней
                    for (i = 3; i<n - 1; i++)
                    {
                        a_inv[i][i - 1] = -J_inv[i - 1];
                        a_inv[i][i] = J_inv[i] + J_inv[i];
                        a_inv[i][i + 1] = -J_inv[i + 1];
                    }
                    //последняя (n-1)
                    if (n>3)
                    {
                        a_inv[n - 1][n - 2] = -J_inv[n - 2];
                        a_inv[n - 1][n - 1] = J_inv[n - 1] + J_inv[n - 1];
                    }

                    //(TJ^{-1})T^{-1}
                    //первая строка
                    aa[1] = a_inv[0][0];  aa[2] = a_inv[0][1];
                    a_inv[0][0] = aa[1] * (double)n + aa[2] * (double)(n - 1);
                    s = aa[1] + aa[2];
                    for (j = 1; j<n; j++) a_inv[0][j] = s*(double)(n - j);
                    //до последней
                    for (i = 1; i<n - 1; i++)
                    {
                        aa[0] = a_inv[i][i - 1];  aa[1] = a_inv[i][i];  aa[2] = a_inv[i][i + 1];
                        for (j = 0; j<i; j++) a_inv[i][j] = aa[0] * (double)(n - i + 1) + aa[1] * (double)(n - i) + aa[2] * (double)(n - i - 1);
                        s = aa[0] + aa[1];
                        a_inv[i][i] = s*(double)(n - i) + aa[2] * (double)(n - i - 1);
                        s += aa[2];
                        for (j = i + 1; j<n; j++) a_inv[i][j] = s*(double)(n - j);
                    }
                    //последняя (n-1)
                    aa[0] = a_inv[n - 1][n - 2];  aa[1] = a_inv[n - 1][n - 1];
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j<n - 1; j++) a_inv[n - 1][j] = s;
                    a_inv[n - 1][n - 1] = aa[0] + aa[1];


                    break;
            }//schema

            break;

    }//variant

    //______________________________________________________________________

    double norm, norm_inv;

    norm = matr_inf_norm(a, n);
    cout << " ||  A  || = " << norm << endl;

    norm_inv = matr_inf_norm(a_inv, n);
    cout << " ||A_inv|| = " << norm_inv << endl;

    double obusl = norm*norm_inv;
    cout << " obusl = " << obusl << endl;

    //невязка генерации
    double **r = new double*[n];
    for (i = 0; i < n; i++)
        r[i] = new double[n];
    matr_mul(a, a_inv, r, n);
    for (i = 0; i < n; i++) r[i][i] -= 1.;

    norm = matr_inf_norm(r, n);

}

void Q_matrix(double **pDouble, int n, int schema) {
    int i, j;
    double  q;

    double curr, next = 1.;
    for (j = 0; j<n - 1; j++)
    {
        curr = next;
        next += 1.;

        q = 1. / sqrt(curr*next);
        for (i = 0; i <= j; i++) pDouble[i][j] = q;
        pDouble[j + 1][j] = -sqrt(curr / next);
        for (i = j + 2; i<n; i++) pDouble[i][j] = 0.;
    }

    q = 1. / sqrt((double)n);
    for (i = 0; i<n; i++) pDouble[i][n - 1] = q;
}

double matr_inf_norm(double **pDouble, int n) {
    int i, j;
    double s, norm = 0.;

    for (i = 0; i < n; i++)
    {
        for (s = 0., j = 0; j < n; j++) s += fabs(pDouble[i][j]);
        if (s > norm) norm = s;
    }

    return norm;
}

void matr_mul(double **a, double **b, double **c, int n)
{
    int i, j, k;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (c[i][j] = 0., k = 0; k<n; k++) c[i][j] += a[i][k] * b[k][j];

}