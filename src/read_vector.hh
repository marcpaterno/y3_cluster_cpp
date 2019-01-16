#ifndef Y3_CLUSTER_CPP_READ_VECTOR_HH
#define Y3_CLUSTER_CPP_READ_VECTOR_HH

// Helper function to read a vector<double> from a file with the given filename.
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

template <class XFORM>
inline std::vector<double>
read_vector(char const* filename, XFORM xform)
{
    std::string fname = std::string("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/data/")
                        + filename;
    std::ifstream file(fname);
	if (!file) {
        std::string errmsg("Failed to open file: ");
        errmsg += fname;
        throw std::runtime_error(errmsg);
    }
    double tmp;
    std::vector<double> res;
    while (file >> tmp)
        res.push_back(xform(tmp));
    return res;
}

inline std::vector<double>
read_vector(char const* filename)
{
    return read_vector(filename, [](double x) { return x; });
}

#endif
