from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray
cimport numpy as np
cdef class Plane(Surface):
    cpdef _intersection(self,Ray A)
    cpdef np.ndarray normal(self,ri)
    cpdef tuple inplane_vectors(self)
    cdef public np.ndarray u
    cdef public np.ndarray v
    cdef public tuple uv_poly
