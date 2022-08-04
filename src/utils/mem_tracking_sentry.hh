#ifndef Y3_CLUSTER_MEM_TRACK_SENTRY_HH
#define Y3_CLUSTER_MEM_TRACK_SENTRY_HH

#include <optional>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <string>

namespace y3_cluster {

  class MemTrackSentry {
    std::ostream* os_;
    std::string label_;
    std::string module_label_;
    int volume_id_ = -1;
    int grid_point_id_ = -1;
    size_t integrand_dev_mem_footprint = 0;
    
  public:
    explicit MemTrackSentry(std::ostream* os,
			    std::string module_label,
			    std::string label,
			    size_t object_mem_size = 0,
			    int volume_id = -1,
			    int grid_point_id = -1)
      : os_(os), module_label_(module_label),label_(label),volume_id_(volume_id), grid_point_id_(grid_point_id), integrand_dev_mem_footprint(object_mem_size)
    {
      print_free_device_mem("start");
    }
    
    void print_free_device_mem(std::string user_string){
      if (os_) {
	size_t free_physmem, total_physmem;
	auto rc = cudaMemGetInfo(&free_physmem, &total_physmem);
	
	if( rc != cudaSuccess){
	  fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(rc));
	}
	
	*os_ << free_physmem << ","
	     << integrand_dev_mem_footprint << ","
	     << module_label_ << ","
	     << label_ << ","
	     << user_string << ","
	     << volume_id_ << ","
	     << grid_point_id_
	     << std::endl;
      }
      
    }
    
    ~MemTrackSentry(){
      print_free_device_mem("end");
    }
  };
}

#endif
