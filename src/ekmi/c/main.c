#include "ekmi.c"
#include "utils/png_utils.h"
#include <time.h> 


int main() {

    // read image
    ImgData img_data = read_png_file("/Users/kx/desktop/1.png");
    png_bytep* __img = get_channel_one(img_data.img, img_data.height, img_data.width);
    int _N = img_data.height;
    int** _img = (int**)init_mat2d(_N);
    for (int i = 0; i < _N; i++) {
        for (int j = 0; j < _N; j++) {
            _img[i][j] = __img[i][j];
        }
    }

    time_t start, stop;
    int _order = 16;
    double _p = 0.5;
    
    start = time(NULL);
    memo_ekmi_rst = get_ekmi(_order, _img, _N, _p);
    stop = time(NULL);
    
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
    printf("\nTime taken: %lds\n", stop-start);

	return 0;
}












