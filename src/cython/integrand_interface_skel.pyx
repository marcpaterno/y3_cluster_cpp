cdef class Integrand{typename}:
    """
    The Y3 Cluster Cosmology integrand for {typename}. Computes cluster abundances and weak
    lensing shear for each (richness, redshift) bin.
    """
    cdef Gamma_T_Integrand[{typename}] *_igrnd
    def __cinit__(self, sample, vector[double] radii,
                        lo_bins, zo_bins):
        cdef DataBlock *ptr = NULL
        cdef size_t raw_ptr = 0
        if isinstance(sample, db.DataBlock):
            raw_ptr = sample._ptr
            ptr = <DataBlock*> raw_ptr
            self._igrnd = new Gamma_T_Integrand[{typename}](deref(ptr),
                                                               radii,
                                                               _into_integration_ranges(lo_bins),
                                                               _into_integration_ranges(zo_bins))
        else:
            raise TypeError('Needs DataBlock argument, got {{}}'.format(type(sample)))

    def __dealloc__(self):
        if self._igrnd is not NULL:
            del self._igrnd

    def with_bins(self, lo_bins, zo_bins):
        cdef Integrand{typename} out = None
        cdef vector[IntegrationRange] clo_bins = _into_integration_ranges(lo_bins)
        cdef vector[IntegrationRange] czo_bins = _into_integration_ranges(zo_bins)
        out._igrnd = deref(self._igrnd).new_with_bins(clo_bins, czo_bins)
        return out

    def centered(self,
                 double scaled_lo,
                 double scaled_lt,
                 double scaled_zt,
                 double scaled_lnM,
                 double mu):
        """
        Computes the centered integrand at a specific point, returning the
        entire data vector at that point.

        NB: each input is the _Cuba_ input to the integrand, so each must be
        in the range [0, 1]. To figure out what values of [lo, lt, zt, lnM]
        those correspond to, consult `gamma_t.hh`.
        """
        if self._igrnd is NULL:
            raise RuntimeError('Integrand is NULL')
        return deref(self._igrnd).centered(scaled_lo,
                                           scaled_lt,
                                           scaled_zt,
                                           scaled_lnM,
                                           mu)

    def miscentered(self,
                    double scaled_lo,
                    double scaled_lc,
                    double scaled_lt,
                    double scaled_zt,
                    double scaled_R,
                    double scaled_lnM,
                    double mu,
                    double scaled_theta):
        """
        Computes the miscentered integrand at a specific point, returning the
        entire data vector at that point.

        NB: each input is the _Cuba_ input to the integrand, so each must be
        in the range [0, 1]. To figure out what actual values
        those correspond to, consult `gamma_t.hh`.
        """
        if self._igrnd is NULL:
            raise RuntimeError('Integrand is NULL')
        return deref(self._igrnd).miscentered(scaled_lo,
                                              scaled_lc,
                                              scaled_lt,
                                              scaled_zt,
                                              scaled_R,
                                              scaled_lnM,
                                              mu,
                                              scaled_theta)

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

        raise ValueError('invalid integration algorithm: {{}}'.format(algo))

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

        raise ValueError('invalid integration algorithm: {{}}'.format(algo))
