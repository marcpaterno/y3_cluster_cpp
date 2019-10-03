#ifndef _Y3_CLUSTER_CPP_DATABLOCK_READER_HH
#define _Y3_CLUSTER_CPP_DATABLOCK_READER_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/datablock_status.h"
#include "cosmosis/datablock/ndarray.hh"
#include <iostream>
#include <stdexcept>
#include <string>

template <typename T>
T
get_datablock(cosmosis::DataBlock& db, const char* section, const char* value)
{
  T a;
  auto retval = db.get_val<T>(section, value, a);
  if (retval != DBS_SUCCESS) {
    throw std::runtime_error(
      "Error (" + std::string(datablock_status_str(retval)) +
      ") when retrieving: " + std::string(section) + ":" + std::string(value));
  }
  return a;
}

// Default version
template <typename T>
T
get_datablock(cosmosis::DataBlock& db,
              const char* section,
              const char* value,
              T default_)
{
  T a;
  auto retval = db.get_val<T>(section, value, a);
  if (retval != DBS_SUCCESS) {
    return default_;
  }
  return a;
}

template <>
inline cosmosis::ndarray<double>
get_datablock(cosmosis::DataBlock& db, const char* section, const char* value)
{
  datablock_type_t typ;
  auto retval = db.get_type(section, value, typ);
  if (retval != DBS_SUCCESS || typ != DBT_DOUBLEND) {
    if (typ != DBT_DOUBLEND)
      retval = DBS_WRONG_VALUE_TYPE;
    throw std::runtime_error(
      "Error (" + std::string(datablock_status_str(retval)) +
      ") when retrieving: " + std::string(section) + ":" + std::string(value));
  }
  return db.view<cosmosis::ndarray<double>>(section, value);
}

// Default version
template <>
inline cosmosis::ndarray<double>
get_datablock(cosmosis::DataBlock& db,
              const char* section,
              const char* value,
              cosmosis::ndarray<double> default_)
{
  datablock_type_t typ;
  auto retval = db.get_type(section, value, typ);
  if (retval != DBS_SUCCESS || typ != DBT_DOUBLEND) {
    return default_;
  }
  return db.view<cosmosis::ndarray<double>>(section, value);
}

#endif // _Y3_CLUSTER_CPP_DATABLOCK_READER_HH
