// A minimal, dependency-free JSON value, parser and writer.
//
// Deliberately small: the build must work offline, so no third-party JSON library is
// pulled in. It supports exactly what the trajectory format needs -- objects, arrays,
// numbers, strings, booleans and null -- and reports parse errors with a byte offset.
#pragma once

#include <cstddef>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ssk::json {

class Value;
using Object = std::map<std::string, Value>;
using Array = std::vector<Value>;

// Enumerators are named Arr/Obj rather than Array/Object so they do not shadow the
// `Array` and `Object` type aliases declared just above.
enum class Type { Null, Bool, Number, String, Arr, Obj };

class Value {
 public:
  Value() = default;
  Value(bool b) : type_(Type::Bool), num_(b ? 1.0 : 0.0) {}
  Value(double d) : type_(Type::Number), num_(d) {}
  Value(int i) : type_(Type::Number), num_(static_cast<double>(i)) {}
  Value(std::size_t i) : type_(Type::Number), num_(static_cast<double>(i)) {}
  Value(const char* s) : type_(Type::String), str_(s) {}
  Value(std::string s) : type_(Type::String), str_(std::move(s)) {}
  Value(Array values) : type_(Type::Arr), arr_(std::move(values)) {}
  Value(Object fields) : type_(Type::Obj), obj_(std::move(fields)) {}

  [[nodiscard]] Type type() const { return type_; }
  [[nodiscard]] bool is_null() const { return type_ == Type::Null; }
  [[nodiscard]] bool is_number() const { return type_ == Type::Number; }
  [[nodiscard]] bool is_array() const { return type_ == Type::Arr; }
  [[nodiscard]] bool is_object() const { return type_ == Type::Obj; }
  [[nodiscard]] bool is_string() const { return type_ == Type::String; }

  [[nodiscard]] double num() const { return num_; }
  [[nodiscard]] bool boolean() const { return num_ != 0.0; }
  [[nodiscard]] const std::string& str() const { return str_; }
  [[nodiscard]] const Array& arr() const { return arr_; }
  [[nodiscard]] const Object& obj() const { return obj_; }
  Array& arr() { return arr_; }
  Object& obj() { return obj_; }

  // Object lookup; returns nullptr when absent or when this is not an object.
  [[nodiscard]] const Value* find(const std::string& key) const;

  void write(std::ostream& out, int indent = 2, int depth = 0) const;
  [[nodiscard]] std::string dump(int indent = 2) const;

 private:
  Type type_ = Type::Null;
  double num_ = 0.0;
  std::string str_;
  Array arr_;
  Object obj_;
};

// Throws std::runtime_error on malformed input.
Value parse(const std::string& text);
Value parse_stream(std::istream& in);
Value parse_file(const std::string& path);

}  // namespace ssk::json
