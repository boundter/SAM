// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OPTION_HANDLERS_HPP_
#define INCLUDE_SAM_OPTION_HANDLERS_HPP_

#include <fstream>
#include <iostream>
#include <vector>

#include <boost/program_options.hpp>

#include <sam/option/argument_base.hpp>

namespace std {

template<typename T>
  std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
  for (auto item : vec) {
    os << item << " ";
  }
  return os;
}

}  // namespace std

namespace sam {

/*!
  *  Parse the values from the command line arguments to the variables given as
  *  Arguments. This also build the interface of the help message. Should
  *  options be given that not exist in the arguments list an error will be
  *  thrown.
  *
  *  @param argc The total number of command line arguments.
  *  @param argv The given command line arguments.
  *  @param arguments The whole list of possible arguments.
  */
void ParseArguments(int argc, char* argv[],
                    std::vector<std::unique_ptr<ArgumentBase>>& arguments);

/*!
 *  Write the arguments to the given file. The arguments will be written on a
 *  single line starting with a #. They are written in the form long_name=value
 *  and are separated by a space. Single values of a vector are separated by a
 *  comma.
 */
void WriteArgumentsToFile(std::vector<std::unique_ptr<ArgumentBase>>& arguments,
                          std::fstream& file);

// Implementation

void ParseArguments(int argc, char* argv[],
                    std::vector<std::unique_ptr<ArgumentBase>>& arguments) {
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()("help,h", "show this help message");
  for (auto it = arguments.begin(); it != arguments.end(); ++it) {
    (*it)->AddArgument(desc);
  }
  boost::program_options::variables_map vmap;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, desc), vmap);
  boost::program_options::notify(vmap);
  if (vmap.count("help")) {
    std::cout << desc;
    exit(0);
  }
  for (auto it = arguments.begin(); it != arguments.end(); ++it) {
    (*it)->ParseArgument(vmap);
  }
}

void WriteArgumentsToFile(std::vector<std::unique_ptr<ArgumentBase>>& arguments,
                          std::fstream& file) {
  file << "#";
  for (auto it = arguments.begin(); it != arguments.end(); ++it) {
    file << " " << (*it)->GetName() << "=" << (*it)->GetValueAsString();
  }
  file << "\n";
  file.flush();
}

}  // namespace sam

#endif  // INCLUDE_SAM_OPTION_HANDLERS_HPP_
