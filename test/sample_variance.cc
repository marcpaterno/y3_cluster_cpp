#include <chrono>
#include <cubacpp/cubacpp.hh>
#include <gperftools/profiler.h>
#include <iostream>

#include "models/omega_z_sdss.hh"
#include "models/sample_variance.hh"

using y3_cluster::SampleVariance_t,
      y3_cluster::IntegrationRange,
      y3_cluster::OMEGA_Z_SDSS,
      y3_cluster::bessel_polynomial_integral;

void
write_matter_power(std::shared_ptr<y3_cluster::Interp2D const>& matter_power_lin)
{
    std::ofstream output_file;
    output_file.open("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/sample_variance_output/matter_power.csv");
    output_file << "z,k_h,p_k,p_k_k_cubed\n";
    std::vector<double> zs{{0.1, 0.2, 0.3}};
    y3_cluster::IntegrationRange lnk_range{std::log(0.0001), std::log(1.0)};
    for (auto z : zs) {
        for (auto i = 0; i < 101; i++) {
            const double k = std::exp(lnk_range.transform(i / 100.0)),
                         p_k = matter_power_lin->eval(k, z),
                         p_k_k_cubed = k*k*k * p_k;
            output_file << z << ","
                        << k << ","
                        << p_k << ","
                        << p_k_k_cubed << "\n";
        }
    }
    output_file.close();
}

void
write_bessel_sums(const y3_cluster::SampleVariance_t& sv)
{
    std::ofstream output_file;
    output_file.open("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/sample_variance_output/bessel_sums.csv");
    output_file << "k_h,sum_of_bessels_z0_z0,sum_of_bessels_z1_z1,sum_of_bessels_z2_z2,"
                << "sum_of_bessels_z0_z1,sum_of_bessels_z0_z2,sum_of_bessels_z1_z2\n";
    y3_cluster::IntegrationRange lnk_range{std::log(0.0001), std::log(0.8)};
    for (auto i = 0; i < 1001; i++) {
        const double k = std::exp(lnk_range.transform(i / 1000.0));
        std::cout << "Summing over bessels for k = " << k << '\n';
        output_file << k << ","
                    << sv.compute_sum_over_bessels(k, 0, 0) << ','
                    << sv.compute_sum_over_bessels(k, 1, 1) << ','
                    << sv.compute_sum_over_bessels(k, 2, 2) << ','
                    << sv.compute_sum_over_bessels(k, 0, 1) << ','
                    << sv.compute_sum_over_bessels(k, 0, 2) << ','
                    << sv.compute_sum_over_bessels(k, 1, 2) << '\n';
    }
    output_file.close();
}

int
main()
{
    cubacores(0, 0);
    const double h = 0.771358152;
    auto const zz = read_vector("z_da_test.txt");
    // da_arr in h inverse Mpc
    auto const da_arr = read_vector("d_a_test.txt");
    auto da_f = std::make_shared<y3_cluster::Interp1D const>(zz, da_arr);

    std::cout << "zlow,zhigh,k,l,int\n";

    std::vector<double> ks{{0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.05}};
    std::vector<std::pair<double, double>> zbins{{{0.1, 0.3}, {0.3, 0.5}, {0.5, 0.7}}};
    for (auto [zlow, zhigh] : zbins) {
        for (auto k : ks) {
            for (auto l = 0u; l < 8; l++) {
                const double rmin = (1 + zlow) * da_f->eval(zlow) * h,
                             rmax = (1 + zhigh) * da_f->eval(zhigh) * h;

                const double int_value = bessel_polynomial_integral(2, l, k, rmin, rmax);
                const double normalized [[maybe_unused]] = int_value * 3 / (rmax*rmax*rmax - rmin*rmin*rmin);

                std::cout << zlow << ','
                          << zhigh << ','
                          << k << ','
                          << l << ','
                          << int_value * 3.0 / (rmax*rmax*rmax - rmin*rmin*rmin)<< '\n';
            }
        }
    }

    std::cout << "DCom({0.1, 0.3}) = {" << 1.1 * da_f->eval(0.1) * h << ", " << 1.3 * da_f->eval(0.3) * h << "}\n";

    OMEGA_Z_SDSS omega_z;
    std::vector<y3_cluster::IntegrationRange> zs{{0.1, 0.3}, {0.3, 0.5}, {0.5, 0.7}};
    SampleVariance_t sv(omega_z, zs, h);
    write_bessel_sums(sv);

    std::shared_ptr<const y3_cluster::Interp2D> matter_power_lin = std::make_shared<const y3_cluster::Interp2D>(
                                                read_vector("matter_power_lin/k_h.txt"),
                                                read_vector("matter_power_lin/z.txt"),
                                                read_vector("matter_power_lin/p_k.txt"));
    write_matter_power(matter_power_lin);

    for (const auto k : ks)
        std::cout << "sum(" << k << ") = " << sv.compute_sum_over_bessels(k, 0, 0) << '\n';

    for (auto l = 0u; l < 10; l++)
        std::cout << "kay(" << l << ", " << 0.1 << ") = " << y3_cluster::survey_mask_kay(omega_z, l, 0.1) << '\n';

    std::cout << "About to integrate...\n";

    ProfilerStart("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/sample_variance_output/dump.txt");

    auto start = std::chrono::high_resolution_clock::now();
    auto matrix = sv.compute();
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = stop - start;

    ProfilerFlush();
    ProfilerStop();

    std::cout << "Done integrating!\n"
              << "Time elapsed: " << diff.count() << "s\n"
              << "\\sigma_{ij}^2 =\n";
    for (const auto& row : matrix) {
        for (const auto sigma : row)
            std::cout << sigma << ' ';
        std::cout << '\n';
    }

    return 0;
}
