// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OPTIONS_HPP_
#define INCLUDE_SAM_OPTIONS_HPP_

#include <iostream>
#include <vector>
#include <string>

#include <boost/program_options.hpp>

namespace std {
template<typename T>
  std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
  for (auto item : vec) {
    os << item << " ";
  }
  return os;
}

std::string to_string(const std::string& value) {
  return value;
}
}

namespace sam {
template<typename T>
std::string ValueAsStringHelper(T value) {
  return std::to_string(value);
}

template<typename T>
std::string ValueAsStringHelper(std::vector<T> value) {
  std::stringstream ss;
  for (auto it = value.begin(); it != value.end(); ++it) {
    ss << (*it);
    if (it + 1 != value.end()) {
      ss << ",";
    }
  }
  return ss.str();
}

/*!
 * \brief Base struct for the definition of a vector holding all the arguments.
 */
struct ArgumentBase {
  std::string short_name;
  std::string long_name;
  std::string name;
  std::string description;

  /*!
   *  @param long_name long name of the command line arguments.
   *  @param short_name short name of the command line arguments.
   *  @param description description of the command line argument.
   */
  ArgumentBase(std::string long_name, std::string short_name,
      std::string description)
  : short_name(short_name), long_name(long_name), description(description) {
    name = long_name + "," + short_name;
  }


  /*!
   *  @param long_name long name of the command line arguments.
   *  @param description description of the command line argument.
   */
  ArgumentBase(std::string long_name, std::string description)
  : long_name(long_name), name(long_name), description(description) {}

  virtual void AddArgument(boost::program_options::options_description& desc) {}

  virtual void ParseArgument(boost::program_options::variables_map& vmap) {}

  std::string GetName() {
    return long_name;
  }

  virtual std::string GetValueAsString() {return std::string();}
};

template<typename T>
/*!
 * \brief Struct to hold the command line arguments.
 */
struct Argument: public ArgumentBase {
  T& value;
  T default_value;

  /*!
   *  @param long_name long name of the command line arguments.
   *  @param description description of the command line argument.
   *  @param value variable holding the argument.
   *  @param default_value default_value of the argument.
   */
  Argument(std::string long_name, std::string description, T& value,
      T default_value)
  : ArgumentBase(long_name, description), value(value),
  default_value(default_value) {}


  /*!
   *  @param long_name long name of the command line arguments.
   *  @param short_name short name of the command line arguments.
   *  @param description description of the command line argument.
   *  @param value variable holding the argument.
   *  @param default_value default_value of the argument.
   */
  Argument(std::string long_name, std::string short_name,
      std::string description, T& value, T default_value)
  : ArgumentBase(long_name, short_name, description), value(value),
  default_value(default_value) {}


  /*!
   *  Adds the arguments to the description.
   */
  void AddArgument(boost::program_options::options_description& desc) {
    desc.add_options()(name.c_str(),
        boost::program_options::value<T>()->default_value(default_value),
        description.c_str());
  }


  /*!
   *  Parse the arguments and save the value to the variable. If no value was
   *  passed it will set the default value.
   */
  void ParseArgument(boost::program_options::variables_map& vmap) {
    if (vmap.count(long_name)) {
      value = vmap[long_name].as<T>();
    }
  }


  std::string GetValueAsString() {
    return ValueAsStringHelper(value);
  }
};


/*

*/
// Specializations
/*
template<typename T> std::string Argument<std::string>
::GetValueAsString() {
  return value;
}

template<typename T> std::string Argument<std::vector<T>>
::GetValueAsString() {
  std::stringstream ss;
  for (auto it = value.begin(); it != value.end(); ++it) {
    ss << (*it);
    if (it + 1 != value.end()) {
      ss << ",";
    }
  }
  return ss.str();
}
*/
}  // namespace sam

#endif  // INCLUDE_SAM_OPTIONS_HPP_
