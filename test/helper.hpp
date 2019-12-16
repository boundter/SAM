// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef TEST_HELPER_HPP_
#define TEST_HELPER_HPP_

#include <vector>
#include <string>

void FillArgv(std::vector<std::string> args, char* argv[]) {
  for (unsigned int i = 0; i < args.size(); ++i) {
    // simply casting the data like
    // argv[i] = const_cast<char*>(args[i].data());
    // leads to freeing of args and the char* are dangling pointers
    argv[i] = new char[args[i].size()+1];
    args[i].copy(argv[i], args[i].size()+1);
    argv[i][args[i].size()] = '\0';
  }
  argv[args.size()+1] = nullptr;
}

#endif  // TEST_HELPER_HPP_
