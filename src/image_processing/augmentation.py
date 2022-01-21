import cv2
import numpy as np
import random
from .preprocessing import norm_char_size

def augment_affine(img):
    h, w = img.shape

    rands = [(random.random() * 3 - 1.5) / 10 for _ in range(3)]
    stretch_coefs = [1, random.random() * 0.2 + 0.8]
    random.shuffle(stretch_coefs)

    src_tri = [
        [0, 0], 
        [w, 0], 
        [0, h]
    ]

    dest_tri = [
        [0, h*rands[0]], 
        [w*stretch_coefs[0], h*rands[1]], 
        [w*rands[2], h*stretch_coefs[1]]
    ]

    src_tri = np.array(src_tri).astype(np.float32)
    dest_tri = np.array(dest_tri).astype(np.float32)
    warp_mat = cv2.getAffineTransform(src_tri, dest_tri)
    res = cv2.warpAffine(img, warp_mat, (h, w))

    return res


def augment_rotate(img):
    center = (img.shape[1]//2, img.shape[0]//2)
    angle = random.random() * 10 - 6
    rot_mat = cv2.getRotationMatrix2D(center, angle, 1)
    res = cv2.warpAffine(img, rot_mat, (img.shape[1], img.shape[0]))
    return res


def augment_morph(img):
    max_coef = int(128 * 0.03)
    coef = random.randint(0, max_coef)
    if coef == 0:
        return img

    se = np.ones((coef, coef))
    if random.random() > 0.5:
        return cv2.erode(img, se)
    return cv2.dilate(img, se)


def augment_image(img):
    h, w = img.shape
    res = np.pad(img, (int(h*0.2), int(w * 0.2)), 'constant', constant_values=(0, 0))
    res = augment_affine(res)
    res = augment_rotate(res)
    res = augment_morph(res)
    res = (norm_char_size(res, 128, 10) > res.max() / 2).astype('int')
    return res