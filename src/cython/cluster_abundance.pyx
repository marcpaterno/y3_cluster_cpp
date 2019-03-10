# distutils: language=c++

from cluster_abundance cimport IntegrationRange as cIR
from cluster_abundance cimport size_t, ostringstream, cubacores, Cuhre, Vegas, Suave
from cython.operator cimport dereference as deref

import cosmosis.datablock as db
import numpy as np


cdef class PyIntegrationRange:
    cdef cIR *_ir
    def __cinit__(self, double start=0, double end=1):
        self._ir = new cIR(start, end)

    @staticmethod
    cdef from_cpp(cIR ir):
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
                raise TypeError('Not a valid IntegrationRange: {}'.format(obj))
            output.push_back(cIR(a, b))
        elif isinstance(obj, PyIntegrationRange):
            output.push_back(cIR(obj.transform(0.0), obj.transform(1.0)))
        else:
            raise TypeError('Not a valid IntegrationRange: {}'.format(obj))

    return output


# Example stolen (and slightly modified) from here:
# https://stackoverflow.com/questions/45133276/passing-c-vector-to-numpy-through-cython-without-copying-and-taking-care-of-me
cdef class VecWrapper:
    cdef vector[double] vec
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]

    cdef set_data(self, const vector[double]& data):
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
    def __cinit__(self, __allow_empty=False):
        if not __allow_empty:
            raise RuntimeError('Cannot directly construct IntegrandResult')
        self._res = NULL
        self.__lo_ir = None
        self.__zo_ir = None

    def __dealloc__(self):
        if self._results:
            del self._results

    @staticmethod
    cdef create(const Gamma_T_Integrated_Bin_Result& result):
        obj = IntegrandResult(__allow_empty=True)
        obj._res = new Gamma_T_Integrated_Bin_Result(result)

        obj.__lo_ir = PyIntegrationRange.from_cpp(result.lo_ir)
        obj.__zo_ir = PyIntegrationRange.from_cpp(result.zo_ir)

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
        return self.__lo_ir

    @property
    def zo_ir(self):
        return self.__zo_ir

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
        return self._res.N

    @property
    def Nw(self):
        return self._res.Nw

    @property
    def Nb(self):
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
            raise IndexError('IntegrandResults out of bounds: {}'.format(n))
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


cdef class Integrand:
    """
    The Y3 Cluster Cosmology integrand. Computes cluster abundances and weak
    lensing shear for each (richness, redshift) bin.
    """
    cdef Gamma_T_Integrand[DefaultModels] *_igrnd
    def __cinit__(self, sample, vector[double] radii,
                        lo_bins, zo_bins):
        cdef DataBlock *ptr = NULL
        cdef size_t raw_ptr = 0
        if isinstance(sample, db.DataBlock):
            raw_ptr = sample._ptr
            ptr = <DataBlock*> raw_ptr
            self._igrnd = new Gamma_T_Integrand[DefaultModels](deref(ptr),
                                                               radii,
                                                               _into_integration_ranges(lo_bins),
                                                               _into_integration_ranges(zo_bins))
        else:
            raise TypeError('Needs DataBlock argument, got {}'.format(type(sample)))

    def __dealloc__(self):
        if self._igrnd is not NULL:
            del self._igrnd

    def integrate_centered(self, epsrel=1e-3, epsabs=1e-12, maxeval=100000000, algo='Cuhre'):
        cdef Cuhre cuhre
        cdef Vegas vegas
        cdef Suave suave

        if self._igrnd is NULL:
            raise RuntimeError('Attempting to integrate NULL integrand')

        cubacores(0, 0)

        cuhre.maxeval = maxeval
        suave.maxeval = maxeval
        vegas.maxeval = maxeval

        if algo == 'Cuhre':
            return IntegrandResults.create(deref(self._igrnd).integrate_centered(cuhre, epsrel, epsabs))
        if algo == 'Vegas':
            return IntegrandResults.create(deref(self._igrnd).integrate_centered(vegas, epsrel, epsabs))
        if algo == 'Suave':
            return IntegrandResults.create(deref(self._igrnd).integrate_centered(suave, epsrel, epsabs))

        raise ValueError('invalid integration algorithm: {}'.format(algo))

    def integrate_miscentered(self, epsrel=1e-3, epsabs=1e-12, maxeval=100000000, algo='Cuhre'):
        cdef Cuhre cuhre
        cdef Vegas vegas
        cdef Suave suave

        if self._igrnd is NULL:
            raise RuntimeError('Attempting to integrate NULL integrand')

        cubacores(0, 0)

        cuhre.maxeval = maxeval
        suave.maxeval = maxeval
        vegas.maxeval = maxeval

        if algo == 'Cuhre':
            return IntegrandResults.create(deref(self._igrnd).integrate_miscentered(cuhre, epsrel, epsabs))
        if algo == 'Vegas':
            return IntegrandResults.create(deref(self._igrnd).integrate_miscentered(vegas, epsrel, epsabs))
        if algo == 'Suave':
            return IntegrandResults.create(deref(self._igrnd).integrate_miscentered(suave, epsrel, epsabs))

        raise ValueError('invalid integration algorithm: {}'.format(algo))
