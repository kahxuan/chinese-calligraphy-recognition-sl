"""
This file contains general image utils for preprocessing
"""

import numpy as np
import cv2


"""
Get boundary of foreground
@param img  :
@param axis : 0 for vertical, 1 for horizontal
@return     : (lo, hi) s.t img[lo:hi] removes the background
"""
def get_boundary(img, axis):
    
    n = img.shape[axis]

    # find background indices
    emp = np.where(np.sum(img, axis=1 - axis) == 0)[0]
    
    if len(emp) == 0: # no background
        return (0, n - 1)

    # find consecutives
    diff = (emp - np.roll(emp, 1))[1:]
    tmp = np.where(diff != 1)[0]
    
    if len(tmp) > 0: # background exists on both lower and upper end
        lo = emp[tmp[0]]
        hi = emp[tmp[-1] + 1]

    elif emp[0] == 0: # background only exists on lower end
        lo = emp[-1]
        hi = n - 1
        
    else: # background only exists on upper end
        lo = 0
        hi = emp[0]
    if emp[0] != 0:
        lo = 0
    if emp[-1] != n - 1:
        hi = n - 1

    return (lo, hi)


"""
Remove padding around foreground, leaving only padding of one pixel wide
@param img  : 
@return     : image with border padding removed
"""
def remove_padding(img):
    left, right = get_boundary(img, 1)
    top, bottom = get_boundary(img, 0)
    return img[top:bottom + 1, left:right + 1]


"""
Remove padding around foreground without affecting the center, 
leaving only padding of one pixel wide
ie. paddings removed on the opposite sides are equal
@param img  : 
@return     : image with border padding removed and center retained
"""
def remove_padding_equal(img):
    h, w = img.shape
    left, right = get_boundary(img, 1)
    top, bottom = get_boundary(img, 0)
    
    width_h = min(w - right - 1, left)
    width_v = min(h - bottom - 1, top)

    return img[width_v:h - width_v, width_h:w - width_h]


"""
Pad image s.t (x_c, y_c) is at the center after padding
@param img  : 
@param x_c  :
@param y_c  :
@return     : image with center (x_c, y_c)
"""
def pad_to_center_coord(img, x_c, y_c):
    h, w = img.shape
    
    pad_width = [
        (int(max(h - y_c * 2, 0)), int(max(y_c * 2 - h, 0))),
        (int(max(w - x_c * 2, 0)), int(max(x_c * 2 - w, 0)))
    ]
    
    img = np.pad(img, pad_width, mode='constant', constant_values=0)
    return img


"""
Get image centroid
@param img  : 
@return     : (x_c, y_c), coordinate of centroid
"""
def get_centroid(img):
    
    h, w = img.shape

    m00 = img.sum()

    y = np.repeat(np.expand_dims(np.array(range(h)), 1), w, 1)
    m01 = np.multiply(y, img).sum()

    x = np.repeat(np.expand_dims(np.array(range(w)), 0), h, 0)
    m10 = np.multiply(x, img).sum()
    
    if m00 == 0:
        x_c, y_c = 0, 0
        
    else:
        x_c = int(m10 // m00)
        y_c = int(m01 // m00)

    return x_c, y_c


"""
Pad image to a square image without affecting the center
@param img  : 
@param n    : height/width after padding, must be larger than the initial height/width
@return     : image of equal height and width, and center retained
"""
def pad_to_square(img, n):
    h, w = img.shape
    top = (n - h) // 2
    bottom = (n - h + 1) // 2
    left = (n - w) // 2
    right = (n - w + 1) // 2
    img_square = cv2.copyMakeBorder(
        img, 
        top, bottom, left, right, 
        cv2.BORDER_CONSTANT, value=0
    )
    return img_square