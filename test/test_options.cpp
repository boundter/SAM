// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <string>
#include <vector>

#include "test/catch.hpp"
#include "sam/options.hpp"

TEST_CASE("seting name of option", "[options]") {
  std::string name = "aa";
  std::string description = "Whatever";
  SECTION("ArgumentBase") {
    sam::ArgumentBase base(name, description);
    CHECK(base.GetName() == name);
  }
  SECTION("Argument") {
    double dummy_value;
    sam::Argument<double> argument(name, description, dummy_value, 0.0);
    CHECK(argument.GetName() == name);
  }
}

TEST_CASE("get value of double as string", "[options]") {
  double value = 1.5;
  sam::Argument<double> argument("aa", "foo", value, 0.);
  std::string value_string = std::to_string(value);
  CHECK(argument.GetValueAsString() == value_string);
}

TEST_CASE("get value of int as string", "[options]") {
  int value = 1;
  sam::Argument<int> argument("aa", "foo", value, 0);
  std::string value_string = std::to_string(value);
  CHECK(argument.GetValueAsString() == value_string);
}

TEST_CASE("get value of string as string", "[options]") {
  std::string value = "Test";
  sam::Argument<std::string> argument("aa", "foo", value, "A");
  CHECK(argument.GetValueAsString() == value);
}

TEST_CASE("get value of int vector as string", "[options]") {
  std::vector<int> value({1, 2, 3});
  sam::Argument<std::vector<int>> argument("aa", "foo", value, 
    std::vector<int>({1}));
    std::string value_string = "1,2,3";
    CHECK(argument.GetValueAsString() == value_string);
}

TEST_CASE("get value of double vector as string", "[options]") {
  std::vector<double> value({1.3, 2.5});
  sam::Argument<std::vector<double>> argument("aa", "foo", value, 
    std::vector<double>({1.}));
  std::string value_string = "1.3,2.5";
  CHECK(argument.GetValueAsString() == value_string);
}
