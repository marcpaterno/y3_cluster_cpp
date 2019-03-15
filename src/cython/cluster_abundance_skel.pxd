# distutils: language=c++

from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "<cstddef>" namespace "std":
    cdef cppclass size_t:
        pass

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        pass

cdef extern from "<sstream>" namespace "std":
    cdef cppclass ostringstream:
        string str() const

cdef extern from "/cosmosis/cosmosis/datablock/datablock.hh" namespace "cosmosis":
    cdef cppclass DataBlock:
        pass

cdef extern from "cubacpp.hh" namespace "cubacpp":
    cdef cppclass Cuhre:
        int maxeval
    cdef cppclass Vegas:
        int maxeval
    cdef cppclass Suave:
        int maxeval

cdef extern from "cubacpp.hh":
    void cubacores(int, int)

cdef extern from "integration_range.hh" namespace "y3_cluster":
    cdef cppclass IntegrationRange:
        IntegrationRange(double, double)
        double jacobian() const
        double transform(double) const

cdef extern from "gamma_t.hh" namespace "y3_cluster":
    cdef cppclass Gamma_T_Integrated_Bin_Result:
        Gamma_T_Integrated_Bin_Result(Gamma_T_Integrated_Bin_Result&)
        IntegrationRange lo_ir, zo_ir
        vector[double] radius
        vector[double] gamma_ts, gamma_t_errors, gamma_t_probs
        double N, N_error, N_prob
        double Nw, Nw_error, Nw_prob
        double Nb, Nb_error, Nb_prob

    cdef cppclass Gamma_T_Integrated_Bin_Result_S:
        long long neval
        int nregions, status
        size_t n_richness, n_redshift
        Gamma_T_Integrated_Bin_Result_S(Gamma_T_Integrated_Bin_Result_S&)
        const Gamma_T_Integrated_Bin_Result& operator[](size_t) const
        size_t size()

    cdef ostream& operator<<(ostream& os, Gamma_T_Integrated_Bin_Result_S)
    cdef ostringstream& operator<<(ostringstream& os, Gamma_T_Integrated_Bin_Result_S)

    cdef cppclass Gamma_T_Integrand[Models]:
        Gamma_T_Integrand(DataBlock& sample,
                          vector[double] radii,
                          vector[IntegrationRange] lo_bins,
                          vector[IntegrationRange] zo_bins) except +
        Gamma_T_Integrand(const Gamma_T_Integrand&) except +

        Gamma_T_Integrand *new_with_bins(vector[IntegrationRange] lo_bins, vector[IntegrationRange] zo_bins) except +

        vector[double] centered(double scaled_lo,
                                double scaled_lt,
                                double scaled_zt,
                                double scaled_lnM,
                                double mu) const
        vector[double] miscentered(double scaled_lo,
                                   double scaled_lc,
                                   double scaled_lt,
                                   double scaled_zt,
                                   double scaled_R,
                                   double scaled_lnM,
                                   double mu,
                                   double scaled_theta) const

        Gamma_T_Integrated_Bin_Result_S integrate_centered[I](I, double epsrel, double epsabs)
        Gamma_T_Integrated_Bin_Result_S integrate_miscentered[I](I, double epsrel, double epsabs)
