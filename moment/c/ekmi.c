#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/png_utils.h"
#include "utils/stirling_number1.h"
#include "utils/stirling_number2.h"
#include <time.h> 


int N, M;
double p = 0.5;

png_bytep *img;
double** memo_c;
double** memo_d;
double** memo_K;
double** memo_tilde_K;
double** memo_ekmi_t;
double** memo_ekmi_rst;
stirling_cache1 sc1 = { 0 };
stirling_cache2 sc2 = { 0 };


// int factorial(int x);
// int central_moment(int* img, int p, int q);

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


double central_moment(int x, int y) {
    double retval;
	double bar_x = geometric_moment(1, 0) / geometric_moment(0, 0);
    double bar_y = geometric_moment(0, 1) / geometric_moment(0, 0);
    int i, j;

    retval = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            retval += pow(i - bar_x, x) * pow(j - bar_y, y) * img[i][j];
        }
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


void free_mat2d(double** mat2d, int size) {
    int i;
    for (i = 0; i < size; i++) {
        free(mat2d[i]);
    }
    free(mat2d);
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
            memo_K[n][x] = memo_K[n][x] / (p * (N2 - m));
        }
    }
    return memo_K[n][x];
}


double tilde_K(int n, int m) {
    if (memo_tilde_K[n][m] == 0) {
        int x, y;
        memo_tilde_K[n][m] = 0;
        for (x = 0; x < N; x++) {
            for (y = 0; y < M; y++) {
                memo_tilde_K[n][m] += img[x][y] * K(n, x) * K(m, y);
            }
        }
    }
    return memo_tilde_K[n][m];
}


double c(int n, int i) {
    if (memo_c[n][i] == 0) {
        int k;
        double num, denom;
        memo_c[n][i] = 0;
        for (k = i; k < n + 1; k++) {
            num = pow(-1, k) * factorial(n) * factorial(N - k);
            denom = pow(p, k) * factorial(N) * factorial(k) * factorial(n - k);
            if (denom == 0) {
                denom = 1; // TODO
            }
            memo_c[n][i] += num / denom * stirling_number1(&sc1, k, i);
        }
    }
    return memo_c[n][i];
}


double d(int i, int s) {
    if (memo_d[i][s] == 0) {
        int m;
        double num, denom;
        memo_d[i][s] = 0;
        for (m = s; m < i + 1; m++) {
            num = pow(-1, s) * factorial(m) * factorial(N) * pow(p, m);
            denom = factorial(m - s) * factorial(N - m) * factorial(m);
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

        // TODO: central_moment is 0
        double bar_x = (
                central_moment(0, 0) * tilde_K(1, 0) - 
                central_moment(1, 0) * tilde_K(0, 0)
            ) / (central_moment(1, 1) * tilde_K(0, 0));
        double bar_y = (
                central_moment(0, 0) * tilde_K(0, 1) - 
                central_moment(1, 0) * tilde_K(0, 0)
            ) / (central_moment(1, 1) * tilde_K(0, 0));

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
                                memo_ekmi_t[n][m] += acc6 * tilde_K(u, v);
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

    double term1, term2, term3, term4, term5, term6, term7, term8, term9, term10;
    double acc1, acc2, acc3, acc4, acc5, acc6;

    double lambda = tilde_K(0, 0);
    double u = 2 * c(2, 2) * c(0, 0) / (c(1, 1) * c(1, 1));
    double v = 2 * c(2, 2) * c(1, 0) * c(1, 0) / (c(0, 0) * c(1, 1) * c(1, 1));
    double theta = 0.5 * atan(
            (u * tilde_K(1, 1) - v * tilde_K(0, 0)) / 
            (tilde_K(2, 0) - tilde_K(0, 2))
        );
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
                        term9 = pow(cos(theta), i + t - s);
                        term10 = pow(sin(theta), j - t + s);
                        acc4 = acc3 * term1 * term3 * term9 * term10;
                        for (r = 0; r < i + j - s - t + 1; r++) {
                            term6 = d(i + j - s - t, r);
                            acc5 = acc4 * term6;
                            for (l = 0; l < s + t + 1; l++) {
                                term7 = d(s + t, l);
                                acc6 = acc5 * term7;
                                memo_ekmi_rst[n][m] += acc6 * ekmi_t(r, l);
                            }
                        }
                    }
                }
            }
        }
    }
    return memo_ekmi_rst[n][m];
}


double** get_ekmi(int order, char* fname) {
        
    // read image
    ImgData img_data;
    img_data = read_png_file(fname);
    img = get_channel_one(img_data.img, img_data.height, img_data.width);
    N = img_data.height;
    // M = img_data.width;
    M = N;

    // initialize memo arrays
    int size = fmin(N, M) * 2 + 1;
    memo_c = init_mat2d(size);
    memo_d = init_mat2d(size);
    memo_K = init_mat2d(size);
    memo_tilde_K = init_mat2d(size);
    memo_ekmi_t = init_mat2d(size);
    memo_ekmi_rst = init_mat2d(order);

    // initialize stirling cache
    if (!stirling_cache_create1(&sc1, size)) {
        fprintf(stderr, "Out of memory\n");
        return memo_ekmi_rst;
    }
    if (!stirling_cache_create2(&sc2, size)) {
        fprintf(stderr, "Out of memory\n");
        return memo_ekmi_rst;
    }

    // extract moment
    for (int i=0; i<order; i++) {
        for (int j=0; j<order; j++) {
            ekmi_rst(i, j);
        }
    }
    
    // clean up
    free(img);
    free_mat2d(memo_c, size);
    free_mat2d(memo_d, size);
    free_mat2d(memo_K, size);
    free_mat2d(memo_tilde_K, size);
    free_mat2d(memo_ekmi_t, size);
    stirling_cache_destroy1(&sc1);
    stirling_cache_destroy2(&sc2);

    // free_mat2d(memo_ekmi_rst, order);
    return memo_ekmi_rst;
}



int main() {

    time_t start, stop;
    int order = 5;
    char *fname  =  "/Users/kx/desktop/1.png";

    start = time(NULL);
    memo_ekmi_rst = get_ekmi(order, fname);
    stop = time(NULL);
    
    printf("\n== ekmi_rst  ==\n");
    for (int i=0; i<order; i++) {
        for (int j=0; j<order; j++) {
            printf("%*.4f ", 9, memo_ekmi_rst[i][j]);
        }
        printf("\n");
    }

    free_mat2d(memo_ekmi_rst, order);
    printf("\nTime taken: %lds\n", stop-start);

	return 0;
}












