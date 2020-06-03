// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OPTION_HANDLERS_HPP_
#define INCLUDE_SAM_OPTION_HANDLERS_HPP_

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <boost/program_options.hpp>
#include <boost/lexical_cast/try_lexical_convert.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/value_semantic.hpp>

#include "./argument_base.hpp"

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

// https://stackoverflow.com/questions/4107087/accepting-negative-doubles-with-boostprogram-options
std::vector< boost::program_options::option> ignore_numbers(
    std::vector<std::string>& args) {
    std::vector< boost::program_options::option> result;
    int pos = 0;
    while (!args.empty()) {
        const auto& arg = args[0];
        double num;
        if (boost::conversion::try_lexical_convert(arg, num)) {
            result.push_back(boost::program_options::option());
              boost::program_options::option& opt = result.back();
            opt.position_key = pos++;
            opt.value.push_back(arg);
            opt.original_tokens.push_back(arg);
            args.erase(args.begin());
        } else {
            break;
        }
    }
    return result;
}

void ParseArguments(int argc, char* argv[],
                    std::vector<std::unique_ptr<ArgumentBase>>& arguments) {
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()("help,h", "show this help message");
  for (auto it = arguments.begin(); it != arguments.end(); ++it) {
    (*it)->AddArgument(desc);
  }
  boost::program_options::variables_map vmap;
  boost::program_options::store(
    boost::program_options::command_line_parser(argc, argv)
    .extra_style_parser(&ignore_numbers)
    .options(desc)
    .run(),
    vmap);
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
