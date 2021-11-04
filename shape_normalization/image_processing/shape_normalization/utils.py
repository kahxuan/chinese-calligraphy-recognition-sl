"""
Utils for shape normalization
"""

import numpy as np
from ..preprocessing import get_centroid


"""
Image resampling by bilinear interpolation
If x_mapped and y_mapped are 1D arrays, they will be repeated to form 2D arrays
@param x_mapped  : array of shape (m, n), OR array of shape (n,), specifies the
                   x-coord of the desired pixel value in the original image
@param y_mapped  : array of shape (m, n), OR array of shape (m,), specifies the
                   y-coord of the desired pixel value in the original image
@param img       :
@return          : interpolated image of shape (m, n)
"""
def interpolate_bilinear(x_mapped, y_mapped, img):
    
    # form 2D arrays
    if len(np.array(x_mapped).shape) == 1 and len(np.array(y_mapped).shape) == 1:
        x_mapped, y_mapped = np.meshgrid(x_mapped, y_mapped)
    
    # add background pixels to deal with edge cases
    h, w = img.shape
    img = np.concatenate([img, np.zeros((h, 1))], axis=1)
    img = np.concatenate([img, np.zeros((1, w + 1))], axis=0)
    
    # get neighbouring coordinates
    x_mapped_floor = np.floor(x_mapped).astype(int).clip(min=0, max=w - 1)
    y_mapped_floor = np.floor(y_mapped).astype(int).clip(min=0, max=h - 1)
    x_mapped_ceil = np.ceil(x_mapped).astype(int).clip(min=0, max=w - 1)
    y_mapped_ceil = np.ceil(y_mapped).astype(int).clip(min=0, max=h - 1)
    
    # get weights of each portion
    w_top_left = (1 - (x_mapped % 1)) * (1 - (y_mapped % 1))
    w_top_right = (x_mapped % 1) * (1 - (y_mapped % 1))
    w_btm_left = (1 - (x_mapped % 1)) * (y_mapped % 1)
    w_btm_right = (x_mapped % 1) * (y_mapped % 1)

    # interpolate - sum([area * weight for each portion])
    interpolated = img[y_mapped_floor, x_mapped_floor] * w_top_left + \
        img[y_mapped_ceil, x_mapped_floor] * w_btm_left + \
        img[y_mapped_floor, x_mapped_ceil] * w_top_right + \
        img[y_mapped_ceil, x_mapped_ceil] * w_btm_right
    
    # handle out of bound pixels
    interpolated[x_mapped < 0] = 0
    interpolated[y_mapped < 0] = 0
    interpolated[x_mapped >= w] = 0
    interpolated[y_mapped >= h] = 0

    return interpolated


"""
Compute the second moment (normal or one-sided)
One-sided second moments are used in binoment normalization
@param img       : 
@param one_sided : if True, moment is split into two parts
@return          : second moment (m02, m20) OR 
                   (m02_minus, m02_plus, m20_minus, m20_plus) if one_sided
"""
def get_second_moment(img, one_sided=False):
    
    h, w = img.shape
    x_c, y_c = get_centroid(img)
    
    tmp_m02 = np.multiply(
        np.square((np.array(range(h)) - y_c)), 
        img.sum(axis=1)
    )

    tmp_m20 = np.multiply(
        np.square((np.array(range(w)) - x_c)), 
        img.sum(axis=0)
    )
    
    if not one_sided:
        m02 = tmp_m02.sum() / max(img.sum(), 1e-10)
        m20 = tmp_m20.sum() / max(img.sum(), 1e-10)
        return m02, m20

    else:
        idx = range(y_c)
        m02_minus = tmp_m02[idx].sum() / max(img[idx, :].sum(), 1e-10)
        idx = range(y_c + 1, len(tmp_m02))
        m02_plus = tmp_m02[idx].sum() / max(img[idx, :].sum(), 1e-10)
        idx = range(x_c)
        m20_minus = tmp_m20[idx].sum() / max(img[:, idx].sum(), 1e-10)
        idx = range(x_c + 1, len(tmp_m20))
        m20_plus = tmp_m20[idx].sum() / max(img[:, idx].sum(), 1e-10)

        return m02_minus, m02_plus, m20_minus, m20_plus


"""
Reset image boundary and crop image
@param img       : 
@param bound_x   : (left, right), x-axis boundary
@param bound_y   : (top, bottom), y-axis boundary
@return          : cropped image
"""
def get_bounded_image(img, bound_x, bound_y):
    
    h, w = img.shape
    
    # shift image, together with its boundary

    bound_x_shifted, bound_y_shifted = bound_x, bound_y
    if bound_x[0] < 0:
        bound_x_shifted = (0, bound_x[1] - bound_x[0])
    if bound_y[0] < 0:
        bound_y_shifted = (0, bound_y[1] - bound_y[0])

    pad_width = [
        [max(-bound_y[0], 0), max(bound_y[1] + 1 - h, 0)], 
        [max(-bound_x[0], 0), max(bound_x[1] + 1 - w, 0)]
    ]

    img = np.pad(img, pad_width, mode='constant', constant_values=0)

    # crop image
    img = img[
        bound_y_shifted[0]:bound_y_shifted[1] + 1, 
        bound_x_shifted[0]:bound_x_shifted[1] + 1
    ]
    
    return img