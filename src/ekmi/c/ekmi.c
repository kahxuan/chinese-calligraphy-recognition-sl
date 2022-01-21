#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ekmi.h"
#include "utils/stirling_number1.h"
#include "utils/stirling_number2.h"


int N, M, order;
double p;
int** img;

double** memo_c;
double** memo_d;

double* memo_w;
double* memo_rho;

double** memo_K;
double** memo_tilde_K;
double** memo_ekmi_t;
double** memo_ekmi_rst;
stirling_cache1 sc1 = { 0 };
stirling_cache2 sc2 = { 0 };


double factorial(int x) {  
	double retval = 1;
	if (x > 0) {
		while (x > 0) {
			retval = retval * x;
			x = x - 1;
		}
	}
	else if (x < 0) {
        retval = 0;
	}
	return retval;
}  


int binom(int a, int b) {
    return factorial(a) / (factorial(b) * factorial(a - b));
}


double geometric_moment(int x, int y) {
    double retval;
    int i, j;

    retval = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            retval += pow(i, x) * pow(j, y) * img[i][j];
        }
    }

    return retval;
}


int stirling_number1_signed(int k, int i) {
    int retval = stirling_number1(&sc1, k, i);
    if ((k + i) % 2 == 1) {
        retval = -retval;
    }
    return retval;
}


double** init_mat2d(int size) {
    int i, j;
    double** retval = (double**)malloc(sizeof(double) * size);
    for (i = 0; i < size; i++) {
        retval[i] = (double*)malloc(sizeof(double) * size);
        for (j = 0; j < size; j++) {
            retval[i][j] = 0.f;
        }
    }
    return retval;
}


double* init_mat1d(int size) {
    int i;
    double* retval = (double*)malloc(sizeof(double) * size);
    for (i = 0; i < size; i++) {
        retval[i] = 0.f;
    }
    return retval;
}


void free_mat2d(double** mat2d, int size) {
    int i;
    for (i = 0; i < size; i++) {
        free(mat2d[i]);
    }
    free(mat2d);
}


double w(int x) {
    if (memo_w[x] == 0) {
        int N2 = N - 1;
        if (x == 0) {
            memo_w[x] = pow(1 - p, N2);
        }
        else {
            memo_w[x] = memo_w[x - 1] * (N2 - (x - 1)) * p / (x * (1 - p));
        }
    }
    return memo_w[x];
}


double rho(int n) {
    if (memo_rho[n] == 0) {
        int N2 = N - 1;
        if (n == 0) {
            memo_rho[n] = 1;
        }
        else {
            memo_rho[n] = ( -1 * ((1 - p) / p) * n / (-N2 + n - 1)) * rho(n - 1);
        }
    }
    return memo_rho[n];
}



double K(int n, int x) {
    if (memo_K[n][x] == 0) {
        
        int N2 = N - 1;

        if (n == 0) {
            memo_K[n][x] = 1;
        }
        else if (n == 1) {
            memo_K[n][x] = (1 - (x / (p * N2)));
        }
        else {
            int m = n - 1;
            memo_K[n][x] = 
                (N2 * p - 2 * m * p + m - x) * K(m, x) - 
                    m * (1 - p) * K(m - 1, x);
            memo_K[n][x] = memo_K[n][x] / (p * (N2 - m)); // 0 if order > N / 2
        }
    }
    return memo_K[n][x];
}


double bar_K(int n, int x) {
    double retval;
    retval = sqrt(w(x) / rho(n)) * K(n, x);
    return retval;
}


double tilde_K(int n, int m) {
    // if (memo_tilde_K[n][m] == 0) {
    //     int x, y;
    //     memo_tilde_K[n][m] = 0;
    //     for (x = 0; x < N; x++) {
    //         for (y = 0; y < M; y++) {
    //             // memo_tilde_K[n][m] += img[x][y] * bar_K(n, x) * bar_K(m, y);
    //             // memo_tilde_K[n][m] += img[x][y] * 
    //             //     sqrt(w(x) / rho(n)) * K(n, x) * 
    //             //     sqrt(w(y) / rho(m)) * K(m, y);
    //             memo_tilde_K[n][m] += img[x][y] * K(n, x) * K(m, y);
    //         }
    //     }
    // }
    // return memo_tilde_K[n][m];

    double res = 0;
    int x, y;
    for (x = 0; x < N; x++) {
        for (y = 0; y < M; y++) {
            res += img[x][y] * K(n, x) * K(m, y);
        }
    }
    return res;
}


double tilde_K_t(int n, int m, double x0, double y0) {
    if (memo_tilde_K[n][m] == 0) {
        int x, y;
        memo_tilde_K[n][m] = 0;
        for (x = 0; x < N; x++) {
            for (y = 0; y < M; y++) {
                memo_tilde_K[n][m] += img[x][y] * 
                    // sqrt(w(x - x0) / rho(n)) * K(n, x) * 
                    // sqrt(w(y - y0) / rho(m)) * K(m, y);
                    sqrt(w(x - x0)) * K(n, x) * 
                    sqrt(w(y - y0)) * K(m, y);
                    // sqrt(1 / rho(n)) * K(n, x) * 
                    // sqrt(1 / rho(m)) * K(m, y);
                // memo_tilde_K[n][m] += img[x][y] *  K(n, x) * K(m, y);
                // memo_tilde_K[n][m] += img[x][y] * 
                //     sqrt(w(x) / rho(n)) * K(n, x) * 
                //     sqrt(w(y) / rho(m)) * K(m, y);
            }
        }
    }
    return memo_tilde_K[n][m];
}



double c(int n, int i) {
    int N2 = N - 1;
    if (memo_c[n][i] == 0) {
        int k;
        double num, denom;
        memo_c[n][i] = 0;
        for (k = i; k < n + 1; k++) {
            num = pow(-1, k) * factorial(n) * factorial(N2 - k);
            denom = pow(p, k) * factorial(N2) * factorial(k) * factorial(n - k);
            if (denom == 0) {
                denom = 1; // TODO
            }
            memo_c[n][i] += num / denom * stirling_number1_signed(k, i);
        }
    }
    return memo_c[n][i];
}


double d(int i, int s) {
    int N2 = N - 1;
    if (memo_d[i][s] == 0) {
        int m;
        double num, denom;
        memo_d[i][s] = 0;
        for (m = s; m < i + 1; m++) {
            num = pow(-1, s) * factorial(m) * factorial(N2) * pow(p, m);
            denom = factorial(m - s) * factorial(N2 - m) * factorial(s);
            if (denom == 0) {
                denom = 1; // TODO
            }
            memo_d[i][s] += num / denom * stirling_number2(&sc2, i, m);
        }
    }
    return memo_d[i][s];
}


double ekmi_t(int n, int m) {

    double term1, term2, term3, term4, term5, term6, term7, term8, term9;
    double acc1, acc2, acc3, acc4, acc5, acc6;

    if (memo_ekmi_t[n][m] == 0) {
        memo_ekmi_t[n][m] = 0;

        double bar_x, bar_y;

        bar_x = (
                c(0, 0) * tilde_K(1, 0) - 
                c(1, 0) * tilde_K(0, 0)
            ) / (c(1, 1) * tilde_K(0, 0));
        bar_y = (
                c(0, 0) * tilde_K(0, 1) - 
                c(1, 0) * tilde_K(0, 0)
            ) / (c(1, 1) * tilde_K(0, 0));

        // bar_x = 0;
        // bar_y = 0;
        // if (img[13][13] == 1) {
        //     bar_y += 3;
        // }

        // bar_x = geometric_moment(1, 0) / geometric_moment(0, 0);
        // bar_y = geometric_moment(0, 1) / geometric_moment(0, 0);
        // bar_x = 0;
        // bar_y = 0;

        // printf("%f %f\n", bar_x, bar_y);


        int i, j, s, t, u, v;

        for (i = 0; i < n + 1; i++) {
            term3 = c(n, i);
            acc1 = term3;
            for (j = 0; j < m + 1; j++) {
                term4 = c(m, j);
                acc2 = acc1 * term4;
                for (s = 0; s < i + 1; s++){
                    term1 = binom(i, s);
                    term8 = pow(bar_x, i - s);
                    acc3 = acc2 * term1 * term8;
                    for (t = 0; t < j + 1; t++) {
                        term2 = binom(j, t);
                        term7 = pow(-1, i + j - s - t);
                        term9 = pow(bar_y, j - t);
                        acc4 = acc3 * term2 * term7 * term9;
                        for (u = 0; u < s + 1; u++) {
                            term5 = d(s, u);
                            acc5 = acc4 * term5;
                            for (v = 0; v < t + 1; v++) {
                                term6 = d(t, v);
                                acc6 = acc5 * term6;
                                // memo_ekmi_t[n][m] += acc6 * tilde_K(u, v);
                                memo_ekmi_t[n][m] += acc6 * tilde_K_t(u, v, bar_x, bar_y);
                            }
                        }
                    }
                }
            }
        }
    }
    return memo_ekmi_t[n][m];
}


double ekmi_rst(int n, int m) {

    memo_ekmi_rst[n][m] = ekmi_t(n, m);
    // memo_ekmi_rst[n][m] = tilde_K(n, m);
    return memo_ekmi_rst[n][m];

    double term1, term2, term3, term4, term5, term6, term7, term8, term9, term10;
    double acc1, acc2, acc3, acc4, acc5, acc6;

    double lambda, theta;

    lambda = tilde_K(0, 0);
    double u = 
        2 * c(2, 2) * c(0, 0) / 
        pow(c(1, 1), 2);
    double v = 
        2 * c(2, 2) * pow(c(1, 0), 2) / 
        (c(0, 0) * pow(c(1, 1), 2));
    theta = 0.5 * atan(
            (u * tilde_K(1, 1) - v * tilde_K(0, 0)) / 
            (tilde_K(2, 0) - tilde_K(0, 2))
        );
    
    // theta = 0;
    // // lambda = geometric_moment(0, 0);
    // lambda = 1;

    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    int i, j, s, t, r, l;

    if (memo_ekmi_rst[n][m] == 0) {
        memo_ekmi_rst[n][m] = 0;

        for (i = 0; i < n + 1; i++) {
            term4 = c(n, i);
            acc1 = term4;
            for (j = 0; j < m + 1; j++) {
                term5 = c(m, j);
                term8 = pow(lambda, -(i + j + 2) / 2);
                acc2 = acc1 * term5 * term8;
                for (s = 0; s < i + 1; s++) {
                    term2 = binom(i, s);
                    acc3 = acc2 * term2;
                    for (t = 0; t < j + 1; t++) {
                        term1 = pow(-1, j - t);
                        term3 = binom(j, t);
                        term9 = pow(cos_theta, i + t - s);
                        term10 = pow(sin_theta, j - t + s);
                        acc4 = acc3 * term1 * term3 * term9 * term10;
                        for (r = 0; r < i + j - s - t + 1; r++) {
                            term6 = d(i + j - s - t, r);
                            acc5 = acc4 * term6;
                            for (l = 0; l < s + t + 1; l++) {
                                term7 = d(s + t, l);
                                acc6 = acc5 * term7;
                                // memo_ekmi_rst[n][m] += acc6 * ekmi_t(r, l);
                                memo_ekmi_rst[n][m] += acc6 * tilde_K(r, l);
                            }
                        }
                    }
                }
            }
        }
    }
    return memo_ekmi_rst[n][m];
}


void init_glob(int _order, int** _img, int _N, double _p) {

    int i, j;

    // set global variables
    order = _order;
    N = M = _N;
    p = _p;
    img = (int**)init_mat2d(N);
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            img[i][j] = _img[i][j];
        }
    }

    // initialize memo arrays
    int memo_size = N * 2 + 1;
    memo_c = init_mat2d(memo_size);
    memo_d = init_mat2d(memo_size);
    memo_w = init_mat1d(memo_size);
    memo_rho = init_mat1d(memo_size);
    memo_K = init_mat2d(memo_size);
    memo_tilde_K = init_mat2d(memo_size);
    memo_ekmi_t = init_mat2d(memo_size);
    memo_ekmi_rst = init_mat2d(order);

    // initialize stirling cache
    if (!stirling_cache_create1(&sc1, memo_size)) {
        fprintf(stderr, "Out of memory\n");
        return;
    }
    if (!stirling_cache_create2(&sc2, memo_size)) {
        fprintf(stderr, "Out of memory\n");
        return;
    }
}
    
void free_glob() {

    int memo_size = N * 2 + 1;

    // clean up
    free_mat2d(memo_c, memo_size);
    free_mat2d(memo_d, memo_size);

    free(memo_w);
    free(memo_rho);

    free_mat2d((double**)img, N);
    free_mat2d(memo_K, memo_size);
    free_mat2d(memo_tilde_K, memo_size);
    free_mat2d(memo_ekmi_t, memo_size);

    stirling_cache_destroy1(&sc1);
    stirling_cache_destroy2(&sc2);
}


double** get_ekmi(int _order, int** _img, int _N, double _p) {

    int i, j;

    init_glob(_order, _img, _N, _p);

    // extract moment
    for (i = 0; i < order; i++) {
        for (j = 0; j < order; j++) {
            ekmi_rst(i, j);
        }
    }
    
    free_glob();

    return memo_ekmi_rst;
}


void free_ekmi_rst() {
    free_mat2d(memo_ekmi_rst, order);
}