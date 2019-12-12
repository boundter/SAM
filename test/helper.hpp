// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef TEST_HELPER_HPP_
#define TEST_HELPER_HPP_

#include <vector>
#include <string>

void FillArgv(std::vector<std::string> args, char* argv[]) {
  for (unsigned int i = 0; i < args.size(); ++i) {
    argv[i] = const_cast<char*>(args[i].data());
  }
  argv[args.size()+1] = NULL;
}

#endif  // TEST_HELPER_HPP_
