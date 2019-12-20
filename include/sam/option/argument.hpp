// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OPTION_ARGUMENT_HPP_
#define INCLUDE_SAM_OPTION_ARGUMENT_HPP_

#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include <sam/option/argument_base.hpp>

namespace std {

std::string to_string(const std::string& value) {
  return value;
}

}  // namespace std

namespace sam {

/*!
 * \brief Struct to hold the command line arguments.
 *
 *  The value of the argument will be saved to a given variable of the
 *  appropriate type. If the argument is not specified on the command line the
 *  value will be filled by some given default.
 */
template<typename T>
struct Argument: public ArgumentBase {
  T& value_;
  T default_value_;

  /*!
   *  @param long_name Long name of the command line arguments.
   *  @param description Description of the command line argument.
   *  @param value Variable holding the argument.
   *  @param default_value Default value of the argument.
   */
  explicit Argument(std::string long_name, std::string description, T& value,
                    T default_value);

  /*!
   *  @param long_name Long name of the command line arguments.
   *  @param short_name Short name of the command line arguments.
   *  @param description Description of the command line argument.
   *  @param value Variable holding the argument.
   *  @param default_value Default value of the argument.
   */
  explicit Argument(std::string long_name, std::string short_name,
                    std::string description, T& value, T default_value);

  /*!
   *  Adds the arguments to the description.
   */
  void AddArgument(boost::program_options::options_description& desc) const;

  /*!
   *  Parse the arguments and save the value to the variable. If no value was
   *  passed it will set the default value.
   */
  void ParseArgument(boost::program_options::variables_map& vmap) const;

  /*!
   *  Convert the value to a string for printing.
   */
  std::string GetValueAsString() const;

 private:
  template<typename M>
  void AddArgumentHelper(M default_value,
                         boost::program_options::options_description& desc)
      const;

  template<typename M>
  void AddArgumentHelper(std::vector<M> default_value,
                         boost::program_options::options_description& desc)
      const;

  template<typename M>
  std::string ValueAsStringHelper(M value) const;

  template<typename M>
  std::string ValueAsStringHelper(std::vector<M> value) const;
};

// Implementation

template<typename T>
Argument<T>::Argument(std::string long_name, std::string description, T& value,
                      T default_value)
    : ArgumentBase(long_name, description), value_(value),
      default_value_(default_value) {}

template<typename T>
Argument<T>::Argument(std::string long_name, std::string short_name,
                      std::string description, T& value, T default_value)
    : ArgumentBase(long_name, short_name, description), value_(value),
      default_value_(default_value) {}

template<typename T>
void Argument<T>::AddArgument(
    boost::program_options::options_description& desc) const {
  // Helper function to differentiate between the handling of values and
  // vectors.
  AddArgumentHelper(default_value_, desc);
}

template<typename T>
template<typename M>
void Argument<T>::AddArgumentHelper(
    M default_value, boost::program_options::options_description& desc) const {
  desc.add_options()(name_.c_str(),
                     boost::program_options::value<M>()->
                       default_value(default_value),
                     description_.c_str());
}

template<typename T>
template<typename M>
void Argument<T>::AddArgumentHelper(std::vector<M> default_value,
    boost::program_options::options_description& desc) const {
  desc.add_options()(name_.c_str(),
                      boost::program_options::value<std::vector<M>>()
                        ->multitoken()-> default_value(default_value),
                      description_.c_str());
}

template<typename T>
void Argument<T>::ParseArgument(
      boost::program_options::variables_map& vmap) const {
  if (vmap.count(long_name_)) {
    value_ = vmap[long_name_].as<T>();
  }
}

template<typename T>
std::string Argument<T>::GetValueAsString() const {
  // Helper function to differentiate between the handling of values and
  // vectors.
  return ValueAsStringHelper(value_);
}

template<typename T>
template<typename M>
std::string Argument<T>::ValueAsStringHelper(M value) const {
  return std::to_string(value);
}

template<typename T>
template<typename M>
std::string Argument<T>::ValueAsStringHelper(std::vector<M> value) const {
  std::stringstream ss;
  for (auto it = value.begin(); it != value.end(); ++it) {
    ss << (*it);
    if (it + 1 != value.end()) {
      ss << ",";
    }
  }
  return ss.str();
}

}  // namespace sam

#endif  // INCLUDE_SAM_OPTION_ARGUMENT_HPP_
