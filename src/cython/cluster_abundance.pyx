# distutils: language=c++

from cluster_abundance cimport IntegrationRange as cIR
from cluster_abundance cimport size_t
from cython.operator cimport dereference as deref

import cosmosis.datablock as db

cdef class PyIntegrationRange:
    cdef cIR *_ir
    def __cinit__(self, double start=0, double end=1):
        self._ir = new cIR(start, end)

    def __dealloc__(self):
        del self._ir

    def transform(self, double x):
        return self._ir.transform(x)

    def jacobian(self):
        return self._ir.jacobian()

cdef vector[IntegrationRange] _into_integration_ranges(lst):
    cdef vector[IntegrationRange] output

    for obj in lst:
        if isinstance(obj, tuple):
            a, b = obj
            if not (isinstance(a, float) and isinstance(b, float)):
                raise TypeError('Not a valid IntegrationRange: {}'.format(obj))
            output.push_back(cIR(a, b))
        elif isinstance(obj, PyIntegrationRange):
            output.push_back(cIR(obj.transform(0.0), obj.transform(1.0)))
        else:
            raise TypeError('Not a valid IntegrationRange: {}'.format(obj))

    return output

cdef class Integrand:
    cdef Gamma_T_Integrand[DefaultModels] *_igrnd
    def __cinit__(self, sample, vector[double] radii,
                        lo_bins, zo_bins):
        cdef DataBlock *ptr = NULL
        cdef size_t raw_ptr = 0
        if isinstance(sample, db.DataBlock):
            raw_ptr = sample._ptr
            ptr = <DataBlock*> raw_ptr
            _igrnd = new Gamma_T_Integrand[DefaultModels](deref(ptr), radii,
                                                          _into_integration_ranges(lo_bins),
                                                          _into_integration_ranges(zo_bins))
        else:
            raise TypeError('Needs DataBlock argument, got {}'.format(type(sample)))

    def __dealloc__(self):
        del self._igrnd
