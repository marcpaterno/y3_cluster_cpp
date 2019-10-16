double integrand(double lk, void* params);
double M_to_R(double M, double Omega_m);
int sigma2_at_R_arr(double R, double* k, double* P, int Nk, double* s2);
int sigma2_at_M_arr(double M,
                    double* k,
                    double* P,
                    int Nk,
                    double om,
                    double* s2);
int G_at_sigma_arr(double* sigma, double* G);
int dndM_sigma2_precomputed(double M,
                            double* sigma2,
                            double* sigma2_top,
                            double* sigma2_bot,
                            double Omega_m,
                            double* dndM);
int dndM_at_M_arr(double M,
                  double* k,
                  double* P,
                  int Nk,
                  double om,
                  double* dndM);
double dndM_at_M(double M, double* k, double* P, int Nk, double om);