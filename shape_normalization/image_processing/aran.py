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


# """
# Pad image to a given aspect ratio
# @param img          : 
# @param aspect_ratio : target aspect_ratio
# @return             : image with specified aspect ratio
# """
# def pad_to_aspect_ratio(img, aspect_ratio):

#     h, w = img.shape

#     # update aspect_ratio to w / h
#     if w > h:
#         aspect_ratio = 1 / aspect_ratio

#     # two possible shape: (h, aspect_ratio * h) and (w / aspect_ratio, w)
    
#     diff_h, diff_w = 0, 0
    
#     if aspect_ratio * h > w:
#         diff_w += int((aspect_ratio * h - w) // 2)
        
#     elif w / aspect_ratio > h:
#         diff_h += int((w / aspect_ratio - h) // 2)

#     pad_width = [(diff_h, diff_h), (diff_w, diff_w)]

#     img = np.pad(img, pad_width, mode='constant', constant_values=0)
#     return img
