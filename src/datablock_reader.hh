#ifndef _Y3_CLUSTER_CPP_DATABLOCK_READER_HH
#define _Y3_CLUSTER_CPP_DATABLOCK_READER_HH

#include <stdexcept>
#include <string>
#include <iostream>
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/datablock_status.h"
#include "/cosmosis/cosmosis/datablock/ndarray.hh"

inline const char *
datablock_status_str(DATABLOCK_STATUS dbs)
{
  switch (dbs) {
    case DBS_SUCCESS:
      return "DBS_SUCCESS";
    case DBS_DATABLOCK_NULL:
      return "DBS_DATABLOCK_NULL";
    case DBS_SECTION_NULL:
      return "DBS_SECTION_NULL";
    case DBS_SECTION_NOT_FOUND:
      return "DBS_SECTION_NOT_FOUND";
    case DBS_NAME_NULL:
      return "DBS_NAME_NULL";
    case DBS_NAME_NOT_FOUND:
      return "DBS_NAME_NOT_FOUND";
    case DBS_NAME_ALREADY_EXISTS:
      return "DBS_NAME_ALREADY_EXISTS";
    case DBS_VALUE_NULL:
      return "DBS_VALUE_NULL";
    case DBS_WRONG_VALUE_TYPE:
      return "DBS_WRONG_VALUE_TYPE";
    case DBS_MEMORY_ALLOC_FAILURE:
      return "DBS_MEMORY_ALLOC_FAILURE";
    case DBS_SIZE_NULL:
      return "DBS_SIZE_NULL";
    case DBS_SIZE_NONPOSITIVE:
      return "DBS_SIZE_NONPOSITIVE";
    case DBS_SIZE_INSUFFICIENT:
      return "DBS_SIZE_INSUFFICIENT";
    case DBS_NDIM_NONPOSITIVE:
      return "DBS_NDIM_NONPOSITIVE";
    case DBS_NDIM_OVERFLOW:
      return "DBS_NDIM_OVERFLOW";
    case DBS_NDIM_MISMATCH:
      return "DBS_NDIM_MISMATCH";
    case DBS_EXTENTS_NULL:
      return "DBS_EXTENTS_NULL";
    case DBS_EXTENTS_MISMATCH:
      return "DBS_EXTENTS_MISMATCH";
    case DBS_LOGIC_ERROR:
      return "DBS_LOGIC_ERROR";
    case DBS_USED_DEFAULT:
      return "DBS_USED_DEFAULT";
  }
  return "INVALID DATA BLOCK STATUS";
}

template<typename T>
T
get_datablock(cosmosis::DataBlock& db, const char * section, const char * value)
{
    T a;
    auto retval = db.get_val<T>(section, value, a);
    if (retval != DBS_SUCCESS) {
        throw std::runtime_error("Error (" + std::string(datablock_status_str(retval)) + ") when retrieving: " + std::string(section) + ":" + std::string(value));
    }
    return a;
}

template<> cosmosis::ndarray<double>
get_datablock(cosmosis::DataBlock& db, const char * section, const char * value);

#endif // _Y3_CLUSTER_CPP_DATABLOCK_READER_HH
