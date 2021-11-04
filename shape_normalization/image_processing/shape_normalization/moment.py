"""
Implementation of moment shape normalization as described in the following paper
- Liu, C. & Sako, H. & Fujisawa (2003). Handwritten Chinese Character 
  Recognition: Alternatives to Nonlinear Normalization.
"""

import numpy as np
from .utils import interpolate_bilinear, get_bounded_image, get_second_moment
from ..preprocessing import get_centroid
from ..aran import get_aran, resize_to_aspect_ratio


"""
Compute new image bounds for scaling
@param img     :
@param alpha   : constant alpha in image scaling
@return        : ((left, right), (top, bottom)), new image bounds
"""
def get_mn_scaling(img, alpha=4):
    
    x_c_init, y_c_init = get_centroid(img)
    
    m02, m20 = get_second_moment(img, one_sided=False)
    delta_y = alpha * m02 ** 0.5
    delta_x = alpha * m20 ** 0.5

    delta_y_half = int(delta_y // 2)
    delta_x_half = int(delta_x // 2)
    bound_x = (x_c_init - delta_x_half, x_c_init + delta_x_half)
    bound_y = (y_c_init - delta_y_half, y_c_init + delta_y_half)
    
    return bound_x, bound_y


"""
Compute the coordinate mapping of moment normalization
@param img     :
@param bound_x : (left, right), x-axis boundary 
@param bound_y : (top, bottom), y-axis boundary
@return        : (x_mapped, y_mapped), coordinate mapping of x and y axis
"""
def get_mn_mapping(img, bound_x, bound_y):
        
    x_c_init, y_c_init = get_centroid(img)
    
    h, w = img.shape
    y_c = h // 2
    x_c = w // 2

    delta_x = bound_x[1] - bound_x[0]
    delta_y = bound_y[1] - bound_y[0]

    # coordinate mapping
    x = np.array(range(bound_x[0], bound_x[1]))
    y = np.array(range(bound_y[0], bound_y[1]))
    x_mapped = (x - x_c_init) * w / delta_x + x_c
    x_mapped = np.clip(x_mapped, a_min=0, a_max=w)
    y_mapped = (y - y_c_init) * h / delta_y + y_c
    y_mapped = np.clip(y_mapped, a_min=0, a_max=h)
    
    return x_mapped, y_mapped


"""
Invert coordinate mapping by interpolating original indices with mapped indices
@param x_mapped      : coordinate mapping of x-axis
@param y_mapped      : coordinate mapping of y-axis
@param shape_bounded : shape of image after resetting bounds
@param shape_init    : shape of original image
@return              : inverted coordinate mapping of x and y axis
"""
def inverse_coord_mapping(x_mapped, y_mapped, shape_bounded, shape_init):
    
    h_init, w_init = shape_init
    h_cropped, w_cropped = shape_bounded
    
    scaled_idx = np.linspace(0, h_cropped - 1, h_init)
    y_mapped_inv = interpolate_bilinear(
        y_mapped, [0],
        np.expand_dims(scaled_idx, axis=0)
    )[0, :]

    scaled_idx = np.linspace(0, w_cropped - 1, w_init)
    x_mapped_inv = interpolate_bilinear(
        x_mapped, [0],
        np.expand_dims(scaled_idx, axis=0)
    )[0, :]
    
    return x_mapped_inv, y_mapped_inv


"""
Normalize image using moment normalization, normalized image is resized with 
ARAN
@param img     : 
@param alpha   : constant alpha in image scaling
@return        : normalized image
"""
def normalize(img, alpha=4):
    
    shape_init = img.shape
    
    bound_x, bound_y = get_mn_scaling(img, alpha=alpha)
    x_mapped, y_mapped = get_mn_mapping(img, bound_x, bound_y)
    img = get_bounded_image(img, bound_x, bound_y)
    
    # invert mapping
    shape_bounded = img.shape
    x_mapped_inv, y_mapped_inv = inverse_coord_mapping(
        x_mapped, y_mapped, 
        shape_bounded, shape_init
    )
    
    img = interpolate_bilinear(x_mapped_inv, y_mapped_inv, img)
    img = (img > img.max() / 2 - (img.max()/ 100)).astype('uint8')
    img = resize_to_aspect_ratio(img, get_aran(img))
    
    return img

