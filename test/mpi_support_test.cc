
#include "catch2/catch.hpp"
#include "utils/mpi_support.hh"

TEST_CASE("works_without_mpi")
{
  auto const [local_rank, local_size, global_rank, global_size, initialized] =
    y3_cluster::get_mpi_info();

  CHECK(local_rank == 0);
  CHECK(local_size == 1);
  CHECK(global_rank == 0);
  CHECK(global_size == 1);
  CHECK(initialized == 0);
}
