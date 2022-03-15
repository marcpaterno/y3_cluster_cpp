#include "mpi_support.hh"

#include "mpi.h"
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace y3_cluster {

  mpi_info get_mpi_info() {
    mpi_info result;
    MPI_Initialized(&(result.mpi_initialized));
    if (result.mpi_initialized == 0) return result;

    MPI_Comm_rank(MPI_COMM_WORLD, &(result.global_rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(result.global_size));


    MPI_Comm local;
    MPI_Info info;
    MPI_Info_create(&info);

    int rc = MPI_Comm_split_type(
        MPI_COMM_WORLD,
        MPI_COMM_TYPE_SHARED,
        result.global_rank,
        info,
        &local);
    if (rc != MPI_SUCCESS) throw std::runtime_error("MPI_Comm_split_type failed\n");
    MPI_Comm_rank(local, &(result.local_rank));
    MPI_Comm_size(local, &(result.local_size));
    MPI_Info_free(&info);
    MPI_Comm_free(&local);
    return result;
  }

  std::ostream& operator<<(std::ostream& os, mpi_info const& info) {
    os << "local: " << info.local_rank << '/' << info.local_size
      << " global: " << info.global_rank << '/' << info.global_size;
    return os;
  }

  void record_timestamp(std::ostream& os, char const* label) {
      os << label
        << '\t'
        << std::setprecision(17)
        << MPI_Wtime()
        << std::endl; // we want the flush.
  }



}
