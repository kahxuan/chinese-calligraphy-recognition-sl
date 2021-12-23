// https://gist.github.com/niw/5963798

#include <stdio.h>
#include <stdlib.h>
#include <png.h>

struct imgData {
	int width, height;
	png_bytep *img;
};

typedef struct imgData ImgData;

ImgData read_png_file(char *filename);
png_bytep* get_channel_one(png_bytep *img_full, int height, int width);

ImgData read_png_file(char *filename) {

  ImgData img_data;
  png_byte color_type;
  png_byte bit_depth;
  png_bytep *row_pointers = NULL;

  FILE *fp = fopen(filename, "rb");

  png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if(!png) abort();

  png_infop info = png_create_info_struct(png);
  if(!info) abort();

  if(setjmp(png_jmpbuf(png))) abort();

  png_init_io(png, fp);

  png_read_info(png, info);

  img_data.width      = png_get_image_width(png, info);
  img_data.height     = png_get_image_height(png, info);
  color_type = png_get_color_type(png, info);
  bit_depth  = png_get_bit_depth(png, info);

  // Read any color_type into 8bit depth, RGBA format.
  // See http://www.libpng.org/pub/png/libpng-manual.txt

  if(bit_depth == 16)
    png_set_strip_16(png);

  if(color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_palette_to_rgb(png);

  // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
  if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
    png_set_expand_gray_1_2_4_to_8(png);

  if(png_get_valid(png, info, PNG_INFO_tRNS))
    png_set_tRNS_to_alpha(png);

  // These color_type don't have an alpha channel then fill it with 0xff.
  if(color_type == PNG_COLOR_TYPE_RGB ||
     color_type == PNG_COLOR_TYPE_GRAY ||
     color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

  if(color_type == PNG_COLOR_TYPE_GRAY ||
     color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
    png_set_gray_to_rgb(png);

  png_read_update_info(png, info);

  if (row_pointers) abort();

  row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * img_data.height);
  for(int y = 0; y < img_data.height; y++) {
    row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
  }

  png_read_image(png, row_pointers);

  fclose(fp);

  png_destroy_read_struct(&png, &info, NULL);

  img_data.img = row_pointers;
  return img_data;
}

png_bytep* get_channel_one(png_bytep *img_full, int height, int width) {
  int i, j;
  png_bytep* img = (png_bytep*)malloc(sizeof(png_bytep) * height);
  for(i = 0; i < height; i++) {
    img[i] = (png_byte*)malloc(sizeof(png_bytep) * width);
  }
  for (i = 0; i < height; i++) {
  	for (j = 0; j < width; j ++) {
  		img[i][j] = img_full[i][j * 4] / 255;
  	}
  }
  return img;
}


