#include "catch2/catch.hpp"
#include "utils/timing_sentry.hh"

#include <sstream>

TEST_CASE("null pointer causes no timing recording")
{
  { // Scope only to control lifetime of the sentry.
    y3_cluster::TimingSentry t(nullptr);
  }
}

TEST_CASE("pointer to ostream causes timing recording")
{
  std::ostringstream os;
  CHECK(os.str().empty());
  { // Scope only to control lifetime of the sentry.
    y3_cluster::TimingSentry t(&os);
  }
  CHECK(! os.str().empty());
}

