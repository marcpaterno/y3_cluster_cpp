#ifndef Y3_CLUSTER_CPP_UTILS_MPI_SUPPORT_HH
#define Y3_CLUSTER_CPP_UTILS_MPI_SUPPORT_HH

#include <iosfwd>

namespace y3_cluster {

  // mpi_info contains information regarding the MPI status of the process that
  // calls it.
  struct mpi_info {
    int local_rank = 0;      // rank on this node
    int local_size = 1;      // number of ranks on this node
    int global_rank = 0;     // rank in whole program
    int global_size = 1;     // number of ranks in whole program
    int mpi_initialized = 0; // 0 == MPI not initialized; 1 == MPI initialized
  };

  // Return the MPI information for this rank. We can't tell if the program was
  // launched using MPI, but such a program will not have initialized, and will
  // return a default-constructed mpi_info object.
  mpi_info get_mpi_info();

  std::ostream& operator<<(std::ostream& os, mpi_info const& info);

  void record_timestamp(std::ostream& os, char const* label);
}
#endif
