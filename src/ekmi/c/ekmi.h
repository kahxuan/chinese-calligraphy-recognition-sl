#ifndef EKMI_H_
#define EKMI_H_

double factorial(int x);
int binom(int a, int b);
// double geometric_moment(int x, int y);
double central_moment(int x, int y);
double** init_mat2d(int size);
void free_mat2d(double** mat2d, int size);
double K(int n, int x);
double tilde_K(int n, int m);
double c(int n, int i);
double ekmi_t(int n, int m);
double ekmi_rst(int n, int m);
double** get_ekmi(int _order, int** _img, int _N, double _p);
void free_ekmi_rst();

#endif