// The JSON every result passes through: what the writer writes the reader must read back to the
// bit, and a number the reader cannot read whole must not load as a shorter one.
#include "io/json.hpp"
#include "test_main.hpp"

#include <bit>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>

namespace {

using ssk::json::Array;
using ssk::json::Object;
using ssk::json::Value;

std::string bits(double d) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.17g (0x%016llx)", d,
                static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(d)));
  return buf;
}

bool round_trips(double d) {
  return std::bit_cast<std::uint64_t>(ssk::json::parse(Value(d).dump()).num()) ==
         std::bit_cast<std::uint64_t>(d);
}

bool rejects(const std::string& text) {
  try {
    (void)ssk::json::parse(text);
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}

TEST("json/a number reads back to the bit, from the smallest subnormal to the largest double") {
  for (const double d : {0.0, -0.0, 5e-324, -5e-324, 1e-310, DBL_MIN, DBL_MAX, -DBL_MAX, 0.1,
                         62.615877, 1e15, 1e15 - 1, -1e15 + 1, 9007199254740992.0, 1e16, 1e300,
                         -1e-300, 123456789012345.6, 1732012345678.0}) {
    CHECK_MSG(round_trips(d), bits(d) + " written as " + Value(d).dump());
  }
}

TEST("json/a million doubles of random bits read back to the bit") {
  std::mt19937_64 rng(59);
  int failed = 0;
  std::string first;
  for (int i = 0; i < 1000000; ++i) {
    const double d = std::bit_cast<double>(rng());
    if (std::isfinite(d) && !round_trips(d)) {
      if (failed++ == 0) {
        first = bits(d) + " written as " + Value(d).dump();
      }
    }
  }
  CHECK_MSG(failed == 0, std::to_string(failed) + " did not, the first " + first);
}

// A timestamp in milliseconds is a whole number of 13 digits, and the shortest form of a round
// one is 1.7e+12.
TEST("json/a whole number under 1e15 is written as one, and NaN or infinity as null") {
  CHECK(Value(1700000000000.0).dump() == "1700000000000");
  CHECK(Value(-42.0).dump() == "-42");
  CHECK(Value(0.0).dump() == "0");
  CHECK(Value(2.5).dump() == "2.5");
  CHECK(Value(std::nan("")).dump() == "null");
  CHECK(Value(HUGE_VAL).dump() == "null");
  CHECK(Value(-HUGE_VAL).dump() == "null");
}

TEST("json/a number with anything after what reads as one is rejected, not shortened") {
  for (const char* text : {"1.2.3", "1-2", "1e", "1e+", "--1", "-", "+", "1.5e3.2", "[1.2.3]",
                           "{\"x\": 4-5}", "1e999", "-1e999"}) {
    CHECK_MSG(rejects(text), std::string(text) + " was accepted");
  }
  for (const char* text : {"0", "-0", "1.5", "-2.25e-3", "1E5", "1e+5", "3e-324"}) {
    CHECK_MSG(!rejects(text), std::string(text) + " was rejected");
  }
}

TEST("json/strings round-trip quotes, backslashes, control characters and UTF-8") {
  const std::string text = "say \"hi\"\\ \t\n\r\x01\x1f mopsi/7/1400 Zürich \xe2\x82\xac";
  CHECK(ssk::json::parse(Value(text).dump()).str() == text);
  CHECK(ssk::json::parse(R"("ü€\/\b\f")").str() == "\xc3\xbc\xe2\x82\xac/\b\f");
}

TEST("json/documents round-trip, nested, empty and indented or not") {
  const Value doc(Object{
      {"dim", Value(2)},
      {"name", Value("a")},
      {"points", Value(Array{Value(Array{Value(1.5), Value(-2.0)}), Value(Array{})})},
      {"empty", Value(Object{})},
      {"flags", Value(Array{Value(true), Value(false), Value()})},
  });
  for (const int indent : {0, 2}) {
    const std::string once = doc.dump(indent);
    CHECK_MSG(ssk::json::parse(once).dump(indent) == once, once);
  }
  const Value back = ssk::json::parse(doc.dump());
  CHECK(back.find("points")->arr()[0].arr()[1].num() == -2.0);
  CHECK(back.find("points")->arr()[1].arr().empty());
  CHECK(back.find("flags")->arr()[0].boolean());
  CHECK(back.find("flags")->arr()[2].is_null());
  CHECK(back.find("missing") == nullptr);
  CHECK(back.find("dim")->find("x") == nullptr);
}

TEST("json/malformed documents are rejected with where") {
  for (const char* text : {"", "   ", "[1, 2", "[1 2]", "[1,]", "{\"a\": 1,}", "{\"a\" 1}", "{a: 1}",
                           "\"open", "\"bad \\q escape\"", "nul", "tru", "[1] 2", "{} {}"}) {
    CHECK_MSG(rejects(text), std::string("'") + text + "' was accepted");
  }
}

}  // namespace
