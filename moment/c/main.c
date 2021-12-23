#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "png_utils.h"
#include "stirling_number1.h"
#include "stirling_number2.h"

int N, M;
double p = 0.5;
png_bytep *img;
double** memo_c;
double** memo_d;
double* memo_rho;
double* memo_w;
double** memo_K;
double** memo_tilde_K;
double** memo_ekmi_t;
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


double rising_factorial(int x, int i) {  
    double retval = 1;
    int j;
    for (j = 0; j < i; j++) {
        retval = retval * (x + j);
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


double* init_mat1d(int size) {

    int i;
    double* retval = (double*)malloc(sizeof(double) * size);
    for (i = 0; i < size; i++) {
        retval[i] = 0.f;
    }
    return retval;
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


double ekmi_t(n, m) {
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
            // printf("%d\n", i);
            for (j = 0; j < m + 1; j++) {
                for (s = 0; s < i + 1; s++){
                    for (t = 0; t < j + 1; t++) {
                        for (u = 0; u < s + 1; u++) {
                            for (v = 0; v < t + 1; v++) {
                                memo_ekmi_t[n][m] += binom(i, s) * binom(j, t) * 
                                    c(n, i) * c(m, j) * d(s, u) * d(t, v) *
                                    pow(-1, i + j - s - t) * 
                                    pow(bar_x, i - s) * pow(bar_y, j - t) *
                                    tilde_K(u, v);
                            }
                        }
                    }
                }
            }
        }
    }
    return memo_ekmi_t[n][m];
}

double ekmi_rst(n, m) {

    double retval = 0;
    double lambda = tilde_K(0, 0);
    double u = 2 * c(2, 2) * c(0, 0) / (c(1, 1) * c(1, 1));
    double v = 2 * c(2, 2) * c(1, 0) * c(1, 0) / (c(0, 0) * c(1, 1) * c(1, 1));
    double theta = 0.5 * atan(
            (u * tilde_K(1, 1) - v * tilde_K(0, 0)) / 
            (tilde_K(2, 0) - tilde_K(0, 2))
        );
    int i, j, s, t, r, l;

    for (i = 0; i < n + 1; i++) {
        // printf("%d\n", i);
        for (j = 0; j < m + 1; j++) {
            for (s = 0; s < i + 1; s++){
                for (t = 0; t < j + 1; t++) {
                    for (r = 0; r < i + j - s - t + 1; r++) {
                        for (l = 0; l < s + t + 1; l++) {
                            retval += pow(-1, j - t) *
                                binom(i, s) * binom(j, t) * 
                                c(n, i) * c(m, j) * 
                                d(i + j - s - t, r) * d(s + t, l) *
                                pow(lambda, -(i + j + 2) / 2) *
                                pow(cos(theta), i + t - s) *
                                pow(sin(theta), j - t + s) *
                                ekmi_t(r, l);
                        }
                    }
                }
            }
        }
    }

    return retval;
}



int main() {

    ImgData img_data;

    img_data = read_png_file("/Users/kx/desktop/1.png");
    img = get_channel_one(img_data.img, img_data.height, img_data.width);
    M = img_data.width;
    N = img_data.height;
  
   
	int size = fmin(N, M) * 2 + 1;
	memo_c = init_mat2d(size);
    memo_d = init_mat2d(size);
    memo_rho = init_mat1d(size);
    memo_w = init_mat1d(size);
    memo_K = init_mat2d(size);
    memo_tilde_K = init_mat2d(size);
    memo_ekmi_t = init_mat2d(size);

    // initialize stirling cache
    if (!stirling_cache_create1(&sc1, size)) {
        fprintf(stderr, "Out of memory\n");
        return 1;
    }
    if (!stirling_cache_create2(&sc2, size)) {
        fprintf(stderr, "Out of memory\n");
        return 1;
    }

    // ==== debugging prints =====
	// print image 
	// for (int i = 0; i < N; i++) {
	// 	for (int j = 0; j < M; j++) {
	// 		printf("%d ", img[i][j]);
	// 	}
	// 	printf("\n");
	// }

    // printf("\n====== K ======\n");
    // for (int i=0; i<10; i++) {
    //     for (int j=0; j<10; j++) {
    //         printf("%*.10f ", 13, K(i, j));
    //     }
    //     printf("\n");
    // }

    printf("\n=== tilde_K ===\n");
    for (int i=0; i<30; i++) {
        for (int j=0; j<30; j++) {
            printf("%*.5f ", 10, tilde_K(i, j));
        }
        printf("\n");
    }

    // printf("\n=== ekmi_t  ===\n");
    // for (int i=0; i<10; i++) {
    //     for (int j=0; j<10; j++) {
    //         printf("%*.5f ", 10, ekmi_t(i, j));
    //     }
    //     printf("\n");
    // }

    // printf("\n== ekmi_rst  ==\n");
    // for (int i=0; i<20; i++) {
    //     for (int j=0; j<20; j++) {
    //         printf("%*.5f ", 10, ekmi_rst(i, j));
    //     }
    //     printf("\n");
    // }
	
	
	free(img);
    free_mat2d(memo_c, size);
    free_mat2d(memo_d, size);
    free(memo_rho);
    free(memo_w);
    free_mat2d(memo_K, size);
    free_mat2d(memo_tilde_K, size);
    stirling_cache_destroy1(&sc1);
    stirling_cache_destroy2(&sc2);

	return 0;
}












