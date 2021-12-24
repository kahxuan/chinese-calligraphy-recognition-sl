import ctypes
import numpy as np


def ekmi(so_path):

    ekmi = ctypes.cdll.LoadLibrary(so_path)
    _get_ekmi = ekmi.get_ekmi
    _get_ekmi.argtypes = [
        ctypes.c_int, 
        np.ctypeslib.ndpointer(dtype=np.uintp, flags="C", ndim=1), 
        ctypes.c_int,
        ctypes.c_double
    ]
    _get_ekmi.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
    
    def f(order, img, p):

        imgpp = (img.__array_interface__['data'][0] + np.arange(img.shape[0])* img.strides[0]).astype(np.uintp) 
        _res = _get_ekmi(order, imgpp, img.shape[0], p)

        # copy array
        res = np.zeros((order, order))
        for i in range(order):
            for j in range(order):
                res[i, j] = _res[i][j]

        ekmi.free_ekmi_rst() # free memory

        return res

    return f