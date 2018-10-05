#ifndef _Y3_CLUSTER_CPP_DATABLOCK_READER_HH
#define _Y3_CLUSTER_CPP_DATABLOCK_READER_HH

#include <stdexcept>
#include <string>
#include <iostream>
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/datablock_status.h"
#include "/cosmosis/cosmosis/datablock/ndarray.hh"

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
