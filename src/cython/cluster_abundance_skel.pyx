# distutils: language=c++

from cluster_abundance cimport IntegrationRange as cIR
from cluster_abundance cimport size_t, ostringstream, cubacores, Cuhre, Vegas, Suave
from cython.operator cimport dereference as deref

import cosmosis.datablock as db
import numpy as np


{warning}


cdef class PyIntegrationRange:
    cdef cIR *_ir
    def __cinit__(self, double start=0, double end=1):
        self._ir = new cIR(start, end)

    @staticmethod
    cdef create(cIR ir):
        return PyIntegrationRange(ir.transform(0.0), ir.transform(1.0))

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
            a, b = float(a), float(b)
            if not (isinstance(a, float) and isinstance(b, float)):
                raise TypeError('Not a valid IntegrationRange: {{}}'.format(obj))
            output.push_back(cIR(a, b))
        elif isinstance(obj, PyIntegrationRange):
            output.push_back(cIR(obj.transform(0.0), obj.transform(1.0)))
        else:
            raise TypeError('Not a valid IntegrationRange: {{}}'.format(obj))

    return output


# Example stolen (and slightly modified) from here:
# https://stackoverflow.com/questions/45133276/passing-c-vector-to-numpy-through-cython-without-copying-and-taking-care-of-me
cdef class VecWrapper:
    cdef vector[double] vec
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]
    cdef object suboffsets

    cpdef set_data(self, const vector[double]& data):
        cdef Py_ssize_t itemsize = sizeof(double)

        self.vec = data
        self.shape[0] = self.vec.size()
        self.strides[0] = itemsize
        self.suboffsets = [0]

    def __getbuffer__(self, Py_buffer *buf, int flags):
        cdef Py_ssize_t itemsize = sizeof(double)

        buf.buf = <char *>&(self.vec[0])
        # Doubles = 'd'
        buf.format = 'd'
        buf.internal = NULL
        buf.itemsize = itemsize
        buf.len = self.vec.size() * itemsize
        buf.ndim = 1
        buf.obj = self
        # Make it readonly
        buf.readonly = 1
        buf.shape = self.shape
        buf.strides = self.strides

    def __releasebuffer__(self, Py_buffer *buf):
        pass


cdef class IntegrandResult:
    cdef Gamma_T_Integrated_Bin_Result *_res
    cdef object _lo_ir
    cdef object _zo_ir
    cdef object __radius
    cdef object __gamma_ts
    cdef object __gamma_t_errors
    cdef object __gamma_t_probs
    def __cinit__(self, __allow_empty=False):
        if not __allow_empty:
            raise RuntimeError('Cannot directly construct IntegrandResult')
        self._res = NULL
        self._lo_ir = None
        self._zo_ir = None

    def __dealloc__(self):
        if self._res:
            del self._res

    @staticmethod
    cdef create(const Gamma_T_Integrated_Bin_Result& result):
        obj = IntegrandResult(__allow_empty=True)
        obj._res = new Gamma_T_Integrated_Bin_Result(result)

        obj._lo_ir = PyIntegrationRange.create(result.lo_ir)
        obj._zo_ir = PyIntegrationRange.create(result.zo_ir)

        obj.__radius = VecWrapper()
        obj.__radius.set_data(obj._res.radius)

        obj.__gamma_ts = VecWrapper()
        obj.__gamma_ts.set_data(obj._res.gamma_ts)

        obj.__gamma_t_errors = VecWrapper()
        obj.__gamma_t_errors.set_data(obj._res.gamma_t_errors)

        obj.__gamma_t_probs = VecWrapper()
        obj.__gamma_t_probs.set_data(obj._res.gamma_t_probs)

        return obj

    @property
    def lo_ir(self):
        return self._lo_ir

    @property
    def zo_ir(self):
        return self._zo_ir

    @property
    def gamma_ts(self):
        return np.asarray(self.__gamma_ts)

    @property
    def gamma_t_errors(self):
        return np.asarray(self.__gamma_t_errors)

    @property
    def gamma_t_probs(self):
        return np.asarray(self.__gamma_t_probs)

    @property
    def N(self):
        if self._res is not NULL:
            return self._res.N

    @property
    def N_error(self):
        if self._res is not NULL:
            return self._res.N_error

    @property
    def Nw(self):
        if self._res is not NULL:
            return self._res.Nw

    @property
    def Nb(self):
        if self._res is not NULL:
            return self._res.Nb


cdef class IntegrandResults:
    cdef Gamma_T_Integrated_Bin_Result_S *_results
    def __cinit__(self, __allow_empty=False):
        if not __allow_empty:
            raise RuntimeError('Cannot directly construct IntegrandResults')
        self._results = NULL

    def __dealloc__(self):
        if self._results:
            del self._results

    def __getitem__(self, long n):
        if self._results == NULL:
            raise RuntimeError('Uninitialized IntegrandResults')
        if (n < 0) or ((<size_t>n) >= self._results.size()):
            raise IndexError('IntegrandResults out of bounds: {{}}'.format(n))
        return IntegrandResult.create(deref(self._results)[n])

    def __iter__(self):
        for i in range(self._results.size()):
            yield self[i]

    def __str__(self):
        cdef ostringstream out
        out << deref(self._results)
        return out.str()

    @staticmethod
    cdef create(Gamma_T_Integrated_Bin_Result_S results):
        obj = IntegrandResults(__allow_empty=True)
        obj._results = new Gamma_T_Integrated_Bin_Result_S(results)
        return obj
