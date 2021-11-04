"""
This file contains ARAN utils
"""
import math
import numpy as np
import cv2
from .preprocessing import remove_padding


"""
Get aspect ratio based on ARAN
@param img          : 
@return             : aspect ratio (0 < aspect ratio < 1)
"""
def get_aran(img):
    img = remove_padding(img)
    shape = sorted(img.shape)
    r_init = shape[0] / shape[1]
    aran = math.sin(math.pi / 2 * r_init) ** 0.5
    return aran


"""
Resize image to a given aspect ratio
@param img          : 
@param aspect_ratio : target aspect_ratio
@return             : image with specified aspect ratio
"""
def resize_to_aspect_ratio(img, aspect_ratio):
    
    h, w = img.shape

    # update aspect_ratio to w / h
    if w > h:
        aspect_ratio = 1 / aspect_ratio

    # two possible shape based on h or w, pick the largest
    shape1 = (h, aspect_ratio * h)
    shape2 = (w / aspect_ratio, w)
    target = shape1 if shape1[0] > shape2[0] else shape2
    
    img = cv2.resize(img, (int(target[1]), int(target[0])))
    
    return img
