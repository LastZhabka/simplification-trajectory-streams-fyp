#include "io/json.hpp"

#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace ssk::json {

const Value* Value::find(const std::string& key) const {
  if (type_ != Type::Obj) {
    return nullptr;
  }
  const auto it = obj_.find(key);
  return (it == obj_.end()) ? nullptr : &it->second;
}

namespace {

void write_string(std::ostream& out, const std::string& s) {
  out << '"';
  for (const char c : s) {
    switch (c) {
      case '"': out << "\\\""; break;
      case '\\': out << "\\\\"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x", c);
          out << buf;
        } else {
          out << c;
        }
    }
  }
  out << '"';
}

void write_number(std::ostream& out, double d) {
  if (!std::isfinite(d)) {
    out << "null";  // JSON has no infinity or NaN
    return;
  }
  if (d == static_cast<double>(static_cast<long long>(d)) && std::fabs(d) < 1e15) {
    out << static_cast<long long>(d);
    return;
  }
  // Shortest representation that reads back as the same double: `%.17g` also round-trips,
  // but writes 62.615876999999998 where 62.615877 is exact and half the bytes -- which a
  // directory full of simplified trajectories notices.
  char buf[40];
  const auto [end, ec] = std::to_chars(buf, buf + sizeof buf, d);
  out << std::string_view(buf, static_cast<std::size_t>(end - buf));
}

void pad(std::ostream& out, int indent, int depth) {
  if (indent > 0) {
    out << '\n' << std::string(static_cast<std::size_t>(indent * depth), ' ');
  }
}

// --- parser ---------------------------------------------------------------------

struct Parser {
  const std::string& s;
  std::size_t i = 0;

  [[noreturn]] void fail(const std::string& msg) const {
    throw std::runtime_error("json: " + msg + " at offset " + std::to_string(i));
  }

  void skip_ws() {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r')) {
      ++i;
    }
  }

  char peek() {
    skip_ws();
    if (i >= s.size()) {
      fail("unexpected end of input");
    }
    return s[i];
  }

  void expect(char c) {
    if (peek() != c) {
      fail(std::string("expected '") + c + "'");
    }
    ++i;
  }

  bool literal(const char* lit) {
    const std::size_t n = std::char_traits<char>::length(lit);
    if (s.compare(i, n, lit) == 0) {
      i += n;
      return true;
    }
    return false;
  }

  std::string parse_string() {
    expect('"');
    std::string out;
    while (true) {
      if (i >= s.size()) {
        fail("unterminated string");
      }
      const char c = s[i++];
      if (c == '"') {
        return out;
      }
      if (c != '\\') {
        out.push_back(c);
        continue;
      }
      if (i >= s.size()) {
        fail("unterminated escape");
      }
      const char e = s[i++];
      switch (e) {
        case '"': out.push_back('"'); break;
        case '\\': out.push_back('\\'); break;
        case '/': out.push_back('/'); break;
        case 'b': out.push_back('\b'); break;
        case 'f': out.push_back('\f'); break;
        case 'n': out.push_back('\n'); break;
        case 'r': out.push_back('\r'); break;
        case 't': out.push_back('\t'); break;
        case 'u': {
          if (i + 4 > s.size()) {
            fail("truncated \\u escape");
          }
          const int cp = std::stoi(s.substr(i, 4), nullptr, 16);
          i += 4;
          // Only the BMP subset the trajectory format needs; encode as UTF-8.
          if (cp < 0x80) {
            out.push_back(static_cast<char>(cp));
          } else if (cp < 0x800) {
            out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
          } else {
            out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
            out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
          }
          break;
        }
        default: fail("unknown escape");
      }
    }
  }

  Value parse_value() {
    const char c = peek();
    if (c == '{') {
      ++i;
      Object o;
      if (peek() == '}') {
        ++i;
        return Value(std::move(o));
      }
      while (true) {
        const std::string key = (peek(), parse_string());
        expect(':');
        o.emplace(key, parse_value());
        const char d = peek();
        if (d == ',') {
          ++i;
          continue;
        }
        if (d == '}') {
          ++i;
          return Value(std::move(o));
        }
        fail("expected ',' or '}'");
      }
    }
    if (c == '[') {
      ++i;
      Array a;
      if (peek() == ']') {
        ++i;
        return Value(std::move(a));
      }
      while (true) {
        a.push_back(parse_value());
        const char d = peek();
        if (d == ',') {
          ++i;
          continue;
        }
        if (d == ']') {
          ++i;
          return Value(std::move(a));
        }
        fail("expected ',' or ']'");
      }
    }
    if (c == '"') {
      return Value(parse_string());
    }
    if (literal("true")) {
      return Value(true);
    }
    if (literal("false")) {
      return Value(false);
    }
    if (literal("null")) {
      return Value();
    }
    // number
    const std::size_t start = i;
    if (i < s.size() && (s[i] == '-' || s[i] == '+')) {
      ++i;
    }
    while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.' ||
                            s[i] == 'e' || s[i] == 'E' || s[i] == '+' || s[i] == '-')) {
      ++i;
    }
    if (start == i) {
      fail("unexpected character");
    }
    try {
      return Value(std::stod(s.substr(start, i - start)));
    } catch (const std::exception&) {
      i = start;
      fail("malformed number");
    }
  }
};

}  // namespace

void Value::write(std::ostream& out, int indent, int depth) const {
  switch (type_) {
    case Type::Null: out << "null"; break;
    case Type::Bool: out << (num_ != 0.0 ? "true" : "false"); break;
    case Type::Number: write_number(out, num_); break;
    case Type::String: write_string(out, str_); break;
    case Type::Arr: {
      if (arr_.empty()) {
        out << "[]";
        break;
      }
      // Arrays of numbers (a point) stay on one line; that keeps trajectory files
      // readable instead of one coordinate per line.
      bool flat = true;
      for (const Value& v : arr_) {
        if (v.type_ == Type::Arr || v.type_ == Type::Obj) {
          flat = false;
          break;
        }
      }
      out << '[';
      for (std::size_t k = 0; k < arr_.size(); ++k) {
        if (k != 0) {
          out << (flat ? ", " : ",");
        }
        if (!flat) {
          pad(out, indent, depth + 1);
        }
        arr_[k].write(out, indent, depth + 1);
      }
      if (!flat) {
        pad(out, indent, depth);
      }
      out << ']';
      break;
    }
    case Type::Obj: {
      if (obj_.empty()) {
        out << "{}";
        break;
      }
      out << '{';
      bool first = true;
      for (const auto& [key, val] : obj_) {
        if (!first) {
          out << ',';
        }
        first = false;
        pad(out, indent, depth + 1);
        write_string(out, key);
        out << ": ";
        val.write(out, indent, depth + 1);
      }
      pad(out, indent, depth);
      out << '}';
      break;
    }
  }
}

std::string Value::dump(int indent) const {
  std::ostringstream ss;
  write(ss, indent, 0);
  return ss.str();
}

Value parse(const std::string& text) {
  Parser p{text};
  Value v = p.parse_value();
  p.skip_ws();
  if (p.i != text.size()) {
    p.fail("trailing characters after value");
  }
  return v;
}

Value parse_stream(std::istream& in) {
  std::ostringstream ss;
  ss << in.rdbuf();
  return parse(ss.str());
}

Value parse_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("json: cannot open '" + path + "'");
  }
  return parse_stream(in);
}

}  // namespace ssk::json
