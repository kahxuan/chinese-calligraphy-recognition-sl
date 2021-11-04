"""
Implementation of pseudo-2D bi-moment shape normalization as described in the 
following paper
- Liu, C. & Marukawa, K. (2005). Pseudo two-dimensional shape normalization 
  methods for handwritten Chinese character recognition.
"""

import numpy as np
from .utils import interpolate_bilinear, get_second_moment
from ..preprocessing import get_centroid
from ..aran import get_aran, resize_to_aspect_ratio


"""
Compute horizontal or vertical weight mapping and weighted images
@param img     :
@param axis    : 0 for horizontal strips, 1 for vertical
@param w0      : constant w0, controls strength of upper and lower strips
@return        : ([f], [weight]), list of weighted images of shape (m, n)
                 and list of weights of shape (m,) or (n,)
"""
def get_weighted_fs(img, axis, w0=0.5):
    
    shape = img.shape
    y_c = shape[axis] // 2

    # weight mapping    

    idx = np.array(range(shape[axis]))
    
    w1 = np.array(range(shape[axis]), dtype='float64')
    w1[idx < y_c] = (y_c - w1[idx < y_c]) / y_c * w0
    w1[idx >= y_c] = 0

    w2 = np.array(range(shape[axis]), dtype='float64')
    w2[idx < y_c] = 1 - (y_c - w2[idx < y_c]) / y_c * w0
    w2[idx >= y_c] = 1 - (w2[idx >= y_c] - y_c) / (shape[axis] - y_c) * w0

    w3 = np.array(range(shape[axis]), dtype='float64')
    w3[idx < y_c] = 0
    w3[idx >= y_c] = (w3[idx >= y_c] - y_c) / (shape[axis] - y_c) * w0

    # truncate weights into [0, 1]
    weights = [np.clip(weight, a_min=0, a_max=1) for weight in [w1, w2, w3]]

    # get weighted images (horizontal or vertical strips)
    fs = []
    for weight in weights:
        weight2d = np.expand_dims(weight, axis=1 - axis) \
            .repeat(shape[1 - axis], axis=1 - axis)
        f = np.multiply(img.astype('float64'), weight2d)
        fs.append(f)
    
    return fs, weights


"""
Compute the coordinate mapping for one strip in pseudo-2D bi-moment normalization
@param f_i        : weighted image, one of the horizontal or vertical strips
@param horizontal : set to True if f_i is a horizontal strip
@param beta       : constant beta in image scaling
@return           : coordinate mapping of the corresponding axis
                    (x-axis if hotizontal, y-axis otherwise)
"""
def get_p2dbmn_submapping(f_i, horizontal, beta):

    axis = 1 if horizontal else 0
    c_init = get_centroid(f_i)[1 - axis]
    n = f_i.shape[axis]

    # compute new image bounds

    moments = get_second_moment(f_i, one_sided=True)
    m02_minus, m02_plus, m20_minus, m20_plus = moments

    if horizontal:
        delta_minus = m20_minus ** 0.5 * beta
        delta_plus = m20_plus ** 0.5 * beta

    else:
        delta_minus = m02_minus ** 0.5 * beta
        delta_plus = m02_plus ** 0.5 * beta

    # curve fitting
    bound = ((c_init - delta_minus), (c_init + delta_plus))
    original = [bound[0], c_init, bound[1]]
    normalized = [0, n // 2, n]
    a, b, c = np.polyfit(normalized, original, 2)

    # evaluate polynomial to get mapping
    coords = np.array(range(n))
    coords_mapped = a * coords ** 2 + b * coords + c
    
    return coords_mapped


"""
Compute the coordinate mapping of pseudo-2D bi-moment normalization
@param img     :
@param beta    : constant beta in image scaling
@param w0      : constant w0, controls strength of upper and lower strips
@return        : (x_mapped, y_mapped), coordinate mapping of x and y axis
"""
def get_p2dbmn_mapping(img, beta, w0):
    
    h, w = img.shape

    # get weights and weighted images
    fs_horizontal, weights_horizontal = get_weighted_fs(img, axis=0, w0=w0)
    fs_vertical, weights_vertical = get_weighted_fs(img, axis=1, w0=w0)

    # get mapping of each strip
    mappings_x = [get_p2dbmn_submapping(f_i, horizontal=True, beta=beta) for f_i in fs_horizontal]
    mappings_y = [get_p2dbmn_submapping(f_i, horizontal=False, beta=beta) for f_i in fs_vertical]

    # combine the three mappings using weighted sum for both horizontal and 
    # vertical strips
    weights_h_exp = [np.expand_dims(weight, axis=1).repeat(w, axis=1) for weight in weights_horizontal]
    mappings_x_exp = [np.expand_dims(mapping, axis=0).repeat(h, axis=0) for mapping in mappings_x]
    x_mapped = \
        weights_h_exp[0] * mappings_x_exp[0] + \
        weights_h_exp[1] * mappings_x_exp[1] + \
        weights_h_exp[2] * mappings_x_exp[2]

    weights_v_exp = [np.expand_dims(weight, axis=0).repeat(h, axis=0) for weight in weights_vertical]
    mappings_y_exp = [np.expand_dims(mapping, axis=1).repeat(w, axis=1) for mapping in mappings_y]
    y_mapped = \
        weights_v_exp[0] * mappings_y_exp[0] + \
        weights_v_exp[1] * mappings_y_exp[1] + \
        weights_v_exp[2] * mappings_y_exp[2]
    
    return x_mapped, y_mapped

    
"""
Normalize image using pseudo-2D bi-moment normalization, normalized image is 
resized with ARAN
@param img     : 
@param beta    : constant beta in image scaling
@param w0      : constant w0, controls strength of upper and lower strips
@return        : normalized image
"""
def normalize(img, beta=2, w0=0.5):
    
    x_mapped, y_mapped = get_p2dbmn_mapping(img, beta=beta, w0=w0)
    
    img  = interpolate_bilinear(x_mapped, y_mapped, img)
    img = (img > img.max() / 2 - (img.max()/ 100)).astype('uint8')
    img = resize_to_aspect_ratio(img, get_aran(img))
        
    return img