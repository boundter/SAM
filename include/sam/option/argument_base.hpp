// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OPTION_ARGUMENT_BASE_HPP_
#define INCLUDE_SAM_OPTION_ARGUMENT_BASE_HPP_

#include <string>

#include <boost/program_options.hpp>


namespace sam {

/*!
 *  \brief Base struct for the definition of a vector holding all the arguments.
 */
struct ArgumentBase {
  std::string short_name_;
  std::string long_name_;
  std::string name_;
  std::string description_;

  /*!
   *  @param long_name Long name of the command line arguments (with --).
   *  @param short_name Short name of the command line arguments (with -).
   *  @param description Description of the command line argument.
   */
  explicit ArgumentBase(std::string long_name, std::string short_name,
               std::string description);

  /*!
   *  @param long_name Long name of the command line arguments (with --).
   *  @param description Description of the command line argument.
   */
  explicit ArgumentBase(std::string long_name, std::string description);

  /*!
   *  Adds the arguments to the description.
   */
  virtual void AddArgument(boost::program_options::options_description& desc)
      const;

  /*!
   *  Parse the arguments and save the value to the variable. If no value was
   *  passed it will set the default value.
   */
  virtual void ParseArgument(boost::program_options::variables_map& vmap) const;

  /*!
   *  Convert the value to a string for printing.
   */
  virtual std::string GetValueAsString() const;

  /*!
   *  Get the long name of the argument.
   */
  std::string GetName() const;
};

// Implementation

ArgumentBase::ArgumentBase(std::string long_name, std::string short_name,
                           std::string description)
    : short_name_(short_name), long_name_(long_name),
      description_(description) {
  name_ = long_name + "," + short_name;
}

ArgumentBase::ArgumentBase(std::string long_name, std::string description)
      : long_name_(long_name), name_(long_name), description_(description) {}

void ArgumentBase::AddArgument(
    boost::program_options::options_description& desc) const {}

void ArgumentBase::ParseArgument(
    boost::program_options::variables_map& vmap) const {}

std::string ArgumentBase::GetValueAsString() const { return std::string(); }

std::string ArgumentBase::GetName() const {
  return long_name_;
}

}  // namespace sam

#endif  // INCLUDE_SAM_OPTION_ARGUMENT_BASE_HPP_
