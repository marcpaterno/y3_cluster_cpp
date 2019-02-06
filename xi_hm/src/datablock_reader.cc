#include <datablock_reader.hh>

template<>
cosmosis::ndarray<double>
get_datablock(cosmosis::DataBlock& db, const char * section, const char * value)
{
    datablock_type_t typ;
    auto retval = db.get_type(section, value, typ);
    if (retval != DBS_SUCCESS || typ != DBT_DOUBLEND) {
        if (typ != DBT_DOUBLEND)
            retval = DBS_WRONG_VALUE_TYPE;
        throw std::runtime_error("Error (" + std::string(datablock_status_str(retval)) + ") when retrieving: " + std::string(section) + ":" + std::string(value));
    }
    return db.view<cosmosis::ndarray<double>>(section, value);
}
