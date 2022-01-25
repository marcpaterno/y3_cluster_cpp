#ifndef Y3_CLUSTER_TIMING_SENTRY_HH
#define Y3_CLUSTER_TIMING_SENTRY_HH

#include <chrono>
#include <optional>
#include <ostream>

namespace y3_cluster {

  // TimingSentry is used to (optionally) record timing samples. Create a
  // TimingSentry at a scope to start it running. When it is destroyed, it
  // writes out its timing result to the specified ofstream, if it is in
  // fact present. If the optional<ofstream> is empty, then nothing is
  // written.
  //
  class TimingSentry {
    std::ostream* os_;
    std::chrono::steady_clock::time_point start_;

  public:
    explicit TimingSentry(std::ostream* os)
      : os_(os), start_(std::chrono::steady_clock::now())
    {}

    ~TimingSentry()
    {
      auto stop = std::chrono::steady_clock::now();
      if (os_) {
        auto t = stop - start_;
        auto d =
          std::chrono::duration_cast<std::chrono::microseconds>(t).count();
        *os_ << d << '\n';
      }
    }
  };
}

#endif
