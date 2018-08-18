#include <chrono>
#include <gperftools/profiler.h>
#include <iostream>

#include "sample_variance.hh"

using y3_cluster::EZ,
      y3_cluster::SampleVariance_t,
      y3_cluster::IntegrationRange;

int
main()
{
    cubacores(0, 0);
    auto const zz = read_vector("z_da_test.txt");
    // da_arr in h inverse Mpc
    auto const da_arr = read_vector("d_a_test.txt");
    auto da_f = std::make_shared<y3_cluster::Interp1D const>(zz, da_arr);

    EZ ez{0.3, 0.7, 0};
    SampleVariance_t sv(ez, {da_f, ez, 0.7}, {{{0.1, 0.3}}});//, {0.3, 0.5}, {0.5, 0.7}}});

    ProfilerStart("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/sample_variance_dump.txt");

    std::cout << "About to integrate...\n";

    auto start = std::chrono::high_resolution_clock::now();
    auto matrix = sv();
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
