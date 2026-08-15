#include "test_main.hpp"

#include <cstring>
#include <exception>
#include <iostream>

namespace ssk::test {

std::vector<Case>& registry() {
  static std::vector<Case> cases;
  return cases;
}

int run_all(const char* filter) {
  int passed = 0;
  std::vector<std::string> failures;

  for (const Case& c : registry()) {
    if (filter != nullptr && c.name.find(filter) == std::string::npos) {
      continue;
    }
    try {
      c.fn();
      ++passed;
      std::cout << "  ok   " << c.name << "\n";
    } catch (const Failure& f) {
      failures.push_back(c.name + ": " + f.what);
      std::cout << "  FAIL " << c.name << "\n         " << f.what << "\n";
    } catch (const std::exception& e) {
      failures.push_back(c.name + ": unexpected exception: " + e.what());
      std::cout << "  FAIL " << c.name << "\n         unexpected exception: " << e.what() << "\n";
    }
  }

  std::cout << "\n" << passed << " passed, " << failures.size() << " failed\n";
  return failures.empty() ? 0 : 1;
}

}  // namespace ssk::test

int main(int argc, char** argv) {
  const char* filter = (argc > 1) ? argv[1] : nullptr;
  return ssk::test::run_all(filter);
}
