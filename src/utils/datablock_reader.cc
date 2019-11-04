#include "datablock_reader.hh"

std::vector<double>
get_vector_double(cosmosis::DataBlock& db, const char* section, const char* name)
{
  std::vector<double> doubles;
  auto retval = db.get_val<std::vector<double>>(section, name, doubles) ;
  if (retval == DBS_SUCCESS) return doubles;
  double x;
  retval = db.get_val<double>(section, name, x);
  if (retval == DBS_SUCCESS) return {x};
  std::vector<int> integers;
  retval = db.get_val<std::vector<int>>(section, name, integers);
  if (retval == DBS_SUCCESS) return std::vector<double>(integers.begin(), integers.end());
  int i;
  retval = db.get_val<int>(section, name, i);
  if (retval == DBS_SUCCESS) {
    std::vector<double> result(1);
    result[0] = i;
    return result;
  }
  throw std::runtime_error(std::string("DataBlock contained no parameer in section ") + section + " named " + name + " that could be coerced to type vector<double>");
}


std::vector<double>
get_vector_double(cosmosis::DataBlock& db, std::string section, std::string name)
{
  return get_vector_double(db, section.c_str(), name.c_str());
}
