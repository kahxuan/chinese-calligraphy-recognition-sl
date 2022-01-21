#include "ekmi.c"
#include "utils/png_utils.h"
#include <time.h> 

double poch(n, k) {
    if (n < 0) {
        double retval;
        n = -n;
        retval = pow(-1, k) * factorial(n) / factorial(n - k);
        return retval;
    }
    return 0;
}

double poch_abs(n, k) {    
    double retval;
    retval = factorial(n) / factorial(n - k);
    return retval;
}




int main() {

    // read image
    ImgData img_data = read_png_file("/Users/kx/desktop/1_small_t.png");
    png_bytep* __img = get_channel_one(img_data.img, img_data.height, img_data.width);
    int _N = img_data.height;
    int** _img = (int**)init_mat2d(_N);
    for (int i = 0; i < _N; i++) {
        for (int j = 0; j < _N; j++) {
            _img[i][j] = __img[i][j];
        }
    }

    // time_t start, stop;
    int _order = 7;
    double _p = 0.5;
    

    // ======================== sanity check for rst ==========================

    // int n = 5;
    // int x = 6;
    // int y = 4;

    // double u = 
    //     2 * central_moment(2, 2) * central_moment(0, 0) / 
    //     (central_moment(1, 1) * central_moment(1, 1));
    // double v = 
    //     2 * central_moment(2, 2) * central_moment(1, 0) * central_moment(1, 0) / 
    //     (central_moment(0, 0) * central_moment(1, 1) * central_moment(1, 1));
    // double theta = 0.5 * atan(
    //         (u * tilde_K(1, 1) - v * tilde_K(0, 0)) / 
    //         (tilde_K(2, 0) - tilde_K(0, 2))
    //     );
    // theta = 0;
    // double a11 = cos(theta);
    // double a22 = cos(theta);
    // double a12 = sin(theta);
    // double a21 = -sin(theta);
    // printf("%f %f %f %f %f\n", theta, a11, a22, a12, a21);
    // double res;
    // int i, s;


    // res = K(n, (int)(a11 * x + a12 * y));
    // printf("%.10f\n", res);


    // res = 0;
    // for (k = 0; k < n + 1; k++) {
    //     res += 
    //         pow(-1, k) * poch_abs(n, k) * poch_abs(x, k) / 
    //         (pow(p, k) * poch_abs(N2, k) * factorial(k));
    // }
    // printf("%.10f\n", res);

    // res = 0;
    // for (k = 0; k < n + 1; k++) {
    //     tmp = pow(-1, k) * poch_abs(n, k) / 
    //         (pow(p, k) * poch_abs(N2, k) * factorial(k));
    //     for (i = 0; i < k + 1; i++) {
    //         res += tmp *
    //             stirling_number1_signed(k, i) * pow(x, i);
    //     }
        
    // }
    // printf("%.10f\n", res);

    // res = 0;
    // for (i = 0; i < n + 1; i++) {
    //     res += c(n, i) * pow((int)(a11 * x + a12 * y), i);
    // }
    // printf("%.10f\n", res);


    // res = 0;
    // for (i = 0; i < n + 1; i++) {
    //     for (s = 0; s < i + 1; s++) {
    //         res += binom(i, s) * c(n, i) * pow(a11, i - s) * pow(a12, s) * pow(x, i - s) * pow(y, s);
    //     }
    // }
    // printf("%.10f\n", res);

    // printf("%.10f\n", tilde_K(n, x));
    // printf("%.10f\n", ekmi_rst(n, x));


    // printf("\n====== tilde_K ======\n");
    // for (int i = 0; i < _order; i++) {
    //     for (int j = 0; j < _order; j++) {
    //         printf("%*.5f ", 9, tilde_K(i, j));
    //     }
    //     printf("\n");
    // }

    // printf("\n====== ekmi_rst ======\n");
    // for (int i = 0; i < _order; i++) {
    //     for (int j = 0; j < _order; j++) {
    //         printf("%*.5f ", 9, ekmi_rst(i, j));
    //     }
    //     printf("\n");
    // }


    // free_glob();

    // ========================================================================

    

    // // ========================= sanity check for t ===========================

    init_glob(_order, _img, _N, _p);

    int n = 5;
    int x = 16;
    int x0 = 14;
    int i, s, u;
    double term1, term2, term3, term4, term5;
    double acc1, acc2, acc3;
    double res;
    // int k; 
    double tmp;
    // int N2 = N - 1;

    res = K(n, x - x0);
    printf("%.10f\n", res);

    // res = 0;
    // for (k = 0; k < n + 1; k++) {
    //     res += 
    //         pow(-1, k) * poch_abs(n, k) * poch_abs(x, k) / 
    //         (pow(p, k) * poch_abs(N2, k) * factorial(k));
    // }
    // printf("%.10f\n", res);

    // res = 0;
    // for (k = 0; k < n + 1; k++) {
    //     tmp = pow(-1, k) * poch_abs(n, k) / 
    //         (pow(p, k) * poch_abs(N2, k) * factorial(k));
    //     for (i = 0; i < k + 1; i++) {
    //         res += tmp *
    //             stirling_number1_signed(k, i) * pow(x, i);
    //     }
        
    // }
    // printf("%.10f\n", res);

    res = 0;
    for (i = 0; i < n + 1; i++) {
        tmp = pow(x - x0, i);
        res += c(n, i) * tmp;
    }
    printf("%.10f\n", res);


    res = 0;
    for (i = 0; i < n + 1; i++) {
        for (s = 0; s < i + 1; s++) {
            res += binom(i, s) * c(n, i) * pow(-1, i - s) * pow(x, s) * pow(x0, i - s);    
        }
    }
    printf("%.10f\n", res);


    res = 0;
    for (i = 0; i < n + 1; i++) {
        term2 = c(n, i);
        acc1 = term2;
        for (s = 0; s < i + 1; s++) {
            term1 = binom(i, s);
            term4 = pow(-1, i - s);
            term5 = pow(x0, i - s);
            acc2 = acc1 * term1 * term4 * term5;
            for (u = 0; u < s + 1; u++) {
                term3 = d(s, u);
                acc3 = acc2 * term3;
                res += acc3 * K(u, x);
            }
        }
    }
    printf("%.10f\n", res);

    // printf("\n====== K ======\n");
    // for (int i = 0; i < _order; i++) {
    //     for (int j = 0; j < _order; j++) {
    //         printf("%*.5f ", 9, K(i, j));
    //     }
    //     printf("\n");
    // }

    // printf("\n====== tilde_K ======\n");
    // for (int i = 0; i < _order; i++) {
    //     for (int j = 0; j < _order; j++) {
    //         printf("%*.5f ", 9, tilde_K(i, j));
    //     }
    //     printf("\n");
    // }

    // printf("\n====== ekmi_t ======\n");
    // for (int i = 0; i < _order; i++) {
    //     for (int j = 0; j < _order; j++) {
    //         printf("%*.5f ", 9, ekmi_t(i, j));
    //     }
    //     printf("\n");
    // }


    free_glob();

    // // ========================================================================

    // start = time(NULL);

    memo_ekmi_rst = get_ekmi(_order, _img, _N, _p);
    // stop = time(NULL);
    
    printf("\n== ekmi_rst  ==\n");
    for (int i = 0; i < _order; i++) {
        for (int j = 0; j < _order; j++) {
            printf("%*.4f ", 9, memo_ekmi_rst[i][j]);
        }
        printf("\n");
    }

    // printf("\n====== K ======\n");
    // for (int i = 0; i < _order; i++) {
    //     for (int j = 0; j < _order; j++) {
    //         printf("%*.10f ", 15, K(i, j));
    //     }
    //     printf("\n");
    // }

    // printf("\n======= rho =======\n");
    // for (int i = 0; i < _order; i++) {
    //     printf("%*.10f\n", 15, rho(i));
    // }

    // printf("\n======== w ========\n");
    // for (int i = 0; i < _order; i++) {
    //     printf("%*.50f\n", 55, w(i));
    // }

    free_ekmi_rst();
    free_mat2d((double**)_img, _N);
    free(__img);
    // printf("\nTime taken: %lds\n", stop-start);

	return 0;
}












