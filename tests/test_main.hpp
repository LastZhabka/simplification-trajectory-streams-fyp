// A ~60-line test harness.
//
// Deliberately dependency-free: the build must work offline, and pulling GoogleTest via
// FetchContent would make an offline configure fail. Registration is by static
// initialisation; `ssk_test_main.cpp` provides main().
#pragma once

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace ssk::test {

struct Case {
  std::string name;
  std::function<void()> fn;
};

std::vector<Case>& registry();
int run_all(const char* filter);

struct Registrar {
  Registrar(const char* name, std::function<void()> fn) { registry().push_back({name, std::move(fn)}); }
};

struct Failure {
  std::string what;
};

[[noreturn]] inline void fail(const std::string& msg) { throw Failure{msg}; }

inline void check(bool cond, const std::string& msg) {
  if (!cond) {
    fail(msg);
  }
}

inline void check_close(double a, double b, double tol, const std::string& msg) {
  if (!(std::fabs(a - b) <= tol)) {
    char buf[256];
    std::snprintf(buf, sizeof buf, " (%.17g vs %.17g, tol %.3g)", a, b, tol);
    fail(msg + buf);
  }
}

}  // namespace ssk::test

#define SSK_CONCAT_(a, b) a##b
#define SSK_CONCAT(a, b) SSK_CONCAT_(a, b)

#define TEST(name)                                                            \
  static void SSK_CONCAT(ssk_test_fn_, __LINE__)();                           \
  static const ::ssk::test::Registrar SSK_CONCAT(ssk_test_reg_, __LINE__)(    \
      name, &SSK_CONCAT(ssk_test_fn_, __LINE__));                             \
  static void SSK_CONCAT(ssk_test_fn_, __LINE__)()

#define CHECK(cond) ::ssk::test::check((cond), std::string(#cond) + "  @" __FILE__ ":" + std::to_string(__LINE__))
#define CHECK_MSG(cond, msg) ::ssk::test::check((cond), std::string(msg) + "  @" __FILE__ ":" + std::to_string(__LINE__))
#define CHECK_CLOSE(a, b, tol) ::ssk::test::check_close((a), (b), (tol), std::string(#a " ~= " #b) + "  @" __FILE__ ":" + std::to_string(__LINE__))
