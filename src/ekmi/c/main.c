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
    ImgData img_data = read_png_file("/Users/kx/desktop/test/2_small.png");
    png_bytep* __img = get_channel_one(img_data.img, img_data.height, img_data.width);
    int _N = img_data.height;
    int** _img = (int**)init_mat2d(_N);
    for (int i = 0; i < _N; i++) {
        for (int j = 0; j < _N; j++) {
            _img[i][j] = __img[i][j];
        }
    }

    int _order = 7;
    double _p = 0.5;
    

    memo_ekmi_rst = get_ekmi(_order, _img, _N, _p);
    
    printf("\n== ekmi_rst  ==\n");
    for (int i = 0; i < _order; i++) {
        for (int j = 0; j < _order; j++) {
            printf("%*.4f ", 9, memo_ekmi_rst[i][j]);
        }
        printf("\n");
    }

    free_ekmi_rst();
    free_mat2d((double**)_img, _N);
    free(__img);

	return 0;
}












