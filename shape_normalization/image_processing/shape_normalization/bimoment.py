"""
Implementation of bi-moment shape normalization as described in the following 
paper
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
@param beta    : constant beta in image scaling
@return        : ((left, right), (top, bottom)), new image bounds
"""
def get_bmn_scaling(img, beta=2):
    
    x_c_init, y_c_init = get_centroid(img)

    moments = get_second_moment(img, one_sided=True)
    m02_minus, m02_plus, m20_minus, m20_plus = moments

    delta_y_minus = m02_minus ** 0.5 * beta
    delta_y_plus = m02_plus ** 0.5 * beta
    delta_x_minus = m20_minus ** 0.5 * beta
    delta_x_plus = m20_plus ** 0.5 * beta


    bound_x = (int(x_c_init - delta_x_minus), int(x_c_init + delta_x_plus))
    bound_y = (int(y_c_init - delta_y_minus), int(y_c_init + delta_y_plus))

    return bound_x, bound_y


"""
Compute the coordinate mapping of bi-moment normalization
@param img     :
@param bound_x : (left, right), x-axis boundary 
@param bound_y : (top, bottom), y-axis boundary
@return        : (x_mapped, y_mapped), coordinate mapping of x and y axis
"""
def get_bmn_mapping(img, bound_x, bound_y):
    
    x_c_init, y_c_init = get_centroid(img)
    h, w = img.shape

    # curve fitting for x-axis
    original_x = [bound_x[0], x_c_init, bound_x[1]]
    original_x = [x - bound_x[0] for x in original_x]
    normalized_x = [0, w//2, w]
    a1, b1, c1 = np.polyfit(normalized_x, original_x, 2)

    # curve fitting for y-axis
    original_y = [bound_y[0], y_c_init, bound_y[1]]
    original_y = [y - bound_y[0] for y in original_y]
    normalized_y = [0, h//2, h]
    a2, b2, c2 = np.polyfit(normalized_y, original_y, 2)

    # evaluate polynomial to get mapping
    x = np.array(range(0, w))
    x_mapped = a1 * x ** 2 + b1 * x + c1
    y = np.array(range(0, h))
    y_mapped = a2 * y ** 2 + b2 * y + c2

    # replace critical curves with clipping
    x_mapped = np.clip(x_mapped, a_min=0, a_max=bound_x[1] - bound_x[0])
    y_mapped = np.clip(y_mapped, a_min=0, a_max=bound_y[1] - bound_y[0])
    
    return x_mapped, y_mapped


"""
Normalize image using bi-moment normalization, normalized image is resized with 
ARAN
@param img     : 
@param beta    : constant beta in image scaling
@return        : normalized image
"""
def normalize(img, beta=2):

    x_c_init, y_c_init = get_centroid(img)
    bound_x, bound_y = get_bmn_scaling(img, beta=beta)
    x_mapped, y_mapped = get_bmn_mapping(img, bound_x, bound_y)
    img = get_bounded_image(img, bound_x, bound_y)
    img = interpolate_bilinear(x_mapped, y_mapped, img)
    img = (img > img.max() / 2 - (img.max()/ 100)).astype('uint8')
    
    # resize based on aran
    img = resize_to_aspect_ratio(img, get_aran(img))
    
    return img