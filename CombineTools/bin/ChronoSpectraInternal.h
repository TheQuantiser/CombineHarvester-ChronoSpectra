#ifndef CombineTools_ChronoSpectraInternal_h
#define CombineTools_ChronoSpectraInternal_h

// ChronoSpectra concept and original implementation: Mohammad Abrar Wadud (2024).
// This header is package-private.  It is deliberately kept beside the executable
// rather than exported as part of the CombineHarvester API.

#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "boost/regex.hpp"
#include "TH1F.h"
#include "RooRealVar.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"

namespace ch {
namespace chronospectra {
namespace detail {

inline std::string trim(std::string value) {
  std::string::size_type first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return std::string();
  std::string::size_type last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

inline std::string quote(std::string const& value) {
  return std::string("'") + value + "'";
}

inline bool is_safe_path_component(std::string const& value) {
  if (value.empty() || value == "." || value == "..") return false;
  for (std::string::const_iterator it = value.begin(); it != value.end(); ++it) {
    unsigned char c = static_cast<unsigned char>(*it);
    if (c < 32 || c == 127 || *it == '/' || *it == '\\') return false;
    if (!(std::isalnum(c) || *it == '_' || *it == '-' || *it == '.')) return false;
  }
  return true;
}

inline void require_safe_path_component(std::string const& value,
                                        std::string const& context) {
  if (!is_safe_path_component(value)) {
    throw std::invalid_argument("ChronoSpectra: unsafe " + context + " path component " +
                                quote(value));
  }
}

inline double parse_finite_double(std::string value, std::string const& context) {
  value = trim(value);
  if (value.empty()) {
    throw std::invalid_argument("ChronoSpectra: empty numeric value for " + context);
  }
  errno = 0;
  char *end = nullptr;
  double result = std::strtod(value.c_str(), &end);
  if (errno == ERANGE || end == value.c_str() || *end != '\0' || !std::isfinite(result)) {
    throw std::invalid_argument("ChronoSpectra: expected a finite number for " +
                                context + ", got " + quote(value));
  }
  return result;
}

inline unsigned parse_unsigned(std::string value, std::string const& context) {
  value = trim(value);
  if (value.empty() || value[0] == '-') {
    throw std::invalid_argument("ChronoSpectra: expected an unsigned integer for " + context);
  }
  errno = 0;
  char *end = nullptr;
  unsigned long long result = std::strtoull(value.c_str(), &end, 10);
  if (errno == ERANGE || end == value.c_str() || *end != '\0' ||
      result > static_cast<unsigned long long>(std::numeric_limits<unsigned>::max())) {
    throw std::invalid_argument("ChronoSpectra: expected an unsigned integer for " +
                                context + ", got " + quote(value));
  }
  return static_cast<unsigned>(result);
}

struct NamedGroup {
  std::string name;
  std::vector<std::string> members;
  std::vector<std::string> resolved;
};

struct NamedGroups {
  std::map<std::string, NamedGroup> groups;
  std::vector<std::string> warnings;
};

inline std::vector<std::string> split_trimmed(std::string const& input, char delimiter) {
  std::vector<std::string> result;
  std::string current;
  for (std::string::const_iterator it = input.begin(); it != input.end(); ++it) {
    if (*it == delimiter) {
      result.push_back(trim(current));
      current.clear();
    } else {
      current.push_back(*it);
    }
  }
  result.push_back(trim(current));
  return result;
}

inline NamedGroups parse_named_groups(std::string const& definition,
                                      std::set<std::string> const& available,
                                      std::string const& kind) {
  NamedGroups result;
  if (trim(definition).empty()) return result;

  std::vector<std::string> group_parts = split_trimmed(definition, ';');
  for (std::vector<std::string>::const_iterator group_it = group_parts.begin();
       group_it != group_parts.end(); ++group_it) {
    if (group_it->empty()) {
      throw std::invalid_argument("ChronoSpectra: empty " + kind + " group");
    }
    std::string::size_type colon = group_it->find(':');
    if (colon == std::string::npos || colon != group_it->rfind(':')) {
      throw std::invalid_argument("ChronoSpectra: " + kind + " group must contain exactly one ':' in " +
                                  quote(*group_it));
    }
    NamedGroup group;
    group.name = trim(group_it->substr(0, colon));
    require_safe_path_component(group.name, kind + " group");
    if (result.groups.count(group.name)) {
      throw std::invalid_argument("ChronoSpectra: duplicate " + kind + " group " + quote(group.name));
    }

    std::vector<std::string> members = split_trimmed(group_it->substr(colon + 1), ',');
    std::set<std::string> member_seen;
    for (std::vector<std::string>::const_iterator member_it = members.begin();
         member_it != members.end(); ++member_it) {
      if (member_it->empty()) {
        throw std::invalid_argument("ChronoSpectra: empty member in " + kind + " group " +
                                    quote(group.name));
      }
      if (!member_seen.insert(*member_it).second) {
        result.warnings.push_back("duplicate " + kind + " group member " + quote(*member_it) +
                                  " normalized in " + quote(group.name));
        continue;
      }
      boost::regex expression;
      try {
        expression.assign(*member_it);
      } catch (boost::regex_error const& error) {
        throw std::invalid_argument("ChronoSpectra: invalid regex " + quote(*member_it) +
                                    " in " + kind + " group " + quote(group.name) + ": " +
                                    error.what());
      }
      bool matched = false;
      for (std::set<std::string>::const_iterator name_it = available.begin();
           name_it != available.end(); ++name_it) {
        if (boost::regex_match(*name_it, expression)) {
          matched = true;
          group.resolved.push_back(*name_it);
        }
      }
      if (!matched) {
        throw std::invalid_argument("ChronoSpectra: " + kind + " group member " +
                                    quote(*member_it) + " matches no available name");
      }
      group.members.push_back(*member_it);
    }
    std::sort(group.resolved.begin(), group.resolved.end());
    group.resolved.erase(std::unique(group.resolved.begin(), group.resolved.end()),
                         group.resolved.end());
    if (group.resolved.empty()) {
      throw std::invalid_argument("ChronoSpectra: " + kind + " group " + quote(group.name) +
                                  " resolves to no names");
    }
    result.groups.insert(std::make_pair(group.name, group));
  }
  return result;
}

class OutputManifest {
 public:
  void add(std::string const& path) {
    std::vector<std::string> components = split_trimmed(path, '/');
    if (components.empty() || components.back().empty()) {
      throw std::invalid_argument("ChronoSpectra: empty output path " + quote(path));
    }
    for (std::vector<std::string>::const_iterator it = components.begin();
         it != components.end(); ++it) {
      require_safe_path_component(*it, "output");
    }
    if (!paths_.insert(path).second) {
      throw std::invalid_argument("ChronoSpectra: duplicate output path " + quote(path));
    }
  }

  bool contains(std::string const& path) const { return paths_.count(path) != 0; }
  std::set<std::string> const& paths() const { return paths_; }

 private:
  std::set<std::string> paths_;
};

inline std::size_t checked_square(std::size_t dimension, std::string const& context) {
  if (dimension != 0 && dimension > std::numeric_limits<std::size_t>::max() / dimension) {
    throw std::overflow_error("ChronoSpectra: " + context + " dimension overflows");
  }
  return dimension * dimension;
}

class OnlineCorrelation {
 public:
  explicit OnlineCorrelation(std::size_t dimension)
      : dimension_(dimension), count_(0), mean_(dimension, 0.0),
        co_moment_(checked_square(dimension, "correlation"), 0.0) {}

  void add(std::vector<double> const& values) {
    if (values.size() != dimension_) {
      throw std::invalid_argument("ChronoSpectra: correlation sample dimension mismatch");
    }
    for (std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it) {
      if (!std::isfinite(*it)) {
        throw std::domain_error("ChronoSpectra: non-finite correlation sample");
      }
    }
    ++count_;
    double inverse_count = 1.0 / static_cast<double>(count_);
    std::vector<double> delta(dimension_);
    for (std::size_t i = 0; i < dimension_; ++i) delta[i] = values[i] - mean_[i];
    for (std::size_t i = 0; i < dimension_; ++i) mean_[i] += delta[i] * inverse_count;
    for (std::size_t i = 0; i < dimension_; ++i) {
      for (std::size_t j = 0; j < dimension_; ++j) {
        co_moment_[i * dimension_ + j] += delta[i] * (values[j] - mean_[j]);
      }
    }
  }

  std::size_t count() const { return count_; }
  std::vector<double> const& mean() const { return mean_; }

  std::vector<double> covariance() const {
    std::vector<double> result(co_moment_.size(), 0.0);
    if (count_ == 0) return result;
    double denominator = static_cast<double>(count_);
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = co_moment_[i] / denominator;
    return result;
  }

  std::vector<double> correlation(double roundoff_epsilon = 1e-12) const {
    std::vector<double> covariance_values = covariance();
    std::vector<double> result(covariance_values.size(), 0.0);
    if (count_ == 0) return result;
    std::vector<double> variances(dimension_, 0.0);
    for (std::size_t i = 0; i < dimension_; ++i) {
      double variance = covariance_values[i * dimension_ + i];
      if (variance > 0.0 && std::isfinite(variance)) variances[i] = variance;
    }
    for (std::size_t i = 0; i < dimension_; ++i) {
      for (std::size_t j = 0; j < dimension_; ++j) {
        if (variances[i] == 0.0 || variances[j] == 0.0) continue;
        double value = covariance_values[i * dimension_ + j] /
                       std::sqrt(variances[i] * variances[j]);
        if (!std::isfinite(value)) {
          throw std::domain_error("ChronoSpectra: non-finite correlation result");
        }
        if (value > 1.0 && value <= 1.0 + roundoff_epsilon) value = 1.0;
        if (value < -1.0 && value >= -1.0 - roundoff_epsilon) value = -1.0;
        if (value > 1.0 || value < -1.0) {
          throw std::domain_error("ChronoSpectra: correlation outside [-1,1]");
        }
        result[i * dimension_ + j] = value;
      }
      if (variances[i] > 0.0) result[i * dimension_ + i] = 1.0;
    }
    return result;
  }

 private:
  std::size_t dimension_;
  std::size_t count_;
  std::vector<double> mean_;
  std::vector<double> co_moment_;
};

class ParameterStateGuard {
 public:
  explicit ParameterStateGuard(ch::CombineHarvester& cmb)
      : cmb_(&cmb), snapshot_(cmb.GetParameters()), restored_(false) {}

  ParameterStateGuard(ParameterStateGuard const&) = delete;
  ParameterStateGuard& operator=(ParameterStateGuard const&) = delete;

  void restore() {
    if (restored_) return;
    for (std::vector<ch::Parameter>::const_iterator saved = snapshot_.begin();
         saved != snapshot_.end(); ++saved) {
      ch::Parameter *current = cmb_->GetParameter(saved->name());
      if (!current) {
        throw std::runtime_error("ChronoSpectra: parameter disappeared during restoration: " +
                                 saved->name());
      }
      current->set_frozen(false);
      current->set_range(saved->range_d(), saved->range_u());
      current->set_err_d(saved->err_d());
      current->set_err_u(saved->err_u());
      current->set_val(saved->val());
      ch::Parameter saved_copy(*saved);
      current->groups() = saved_copy.groups();
      current->set_frozen(saved->frozen());
    }
    restored_ = true;
  }

  ~ParameterStateGuard() noexcept {
    if (!restored_) {
      try {
        restore();
      } catch (...) {
      }
    }
  }

 private:
  ch::CombineHarvester *cmb_;
  std::vector<ch::Parameter> snapshot_;
  bool restored_;
};

class ScopedAddDirectory {
 public:
  explicit ScopedAddDirectory(bool enabled)
      : previous_(TH1::AddDirectoryStatus()) {
    TH1::AddDirectory(enabled);
  }
  ~ScopedAddDirectory() { TH1::AddDirectory(previous_); }

 private:
  bool previous_;
};

inline bool compatible_binning(TH1F const& first, TH1F const& second,
                               std::string const& first_name,
                               std::string const& second_name,
                               std::string *reason = nullptr) {
  if (first.GetNbinsX() != second.GetNbinsX()) {
    if (reason) *reason = "regular-bin count differs";
    return false;
  }
  for (int i = 1; i <= first.GetNbinsX(); ++i) {
    std::string left = first.GetXaxis()->GetBinLabel(i);
    std::string right = second.GetXaxis()->GetBinLabel(i);
    if (left != right) {
      if (reason) *reason = "label differs at regular bin " + std::to_string(i) +
                            " (" + quote(left) + " vs " + quote(right) + ")";
      return false;
    }
    double left_edge = first.GetXaxis()->GetBinLowEdge(i);
    double right_edge = second.GetXaxis()->GetBinLowEdge(i);
    double tolerance = 1e-12 + 1e-12 * std::max(std::fabs(left_edge), std::fabs(right_edge));
    if (std::fabs(left_edge - right_edge) > tolerance) {
      if (reason) *reason = "edge differs at regular bin " + std::to_string(i);
      return false;
    }
  }
  double left_edge = first.GetXaxis()->GetBinUpEdge(first.GetNbinsX());
  double right_edge = second.GetXaxis()->GetBinUpEdge(second.GetNbinsX());
  double tolerance = 1e-12 + 1e-12 * std::max(std::fabs(left_edge), std::fabs(right_edge));
  if (std::fabs(left_edge - right_edge) > tolerance) {
    if (reason) *reason = "upper edge differs";
    return false;
  }
  (void)first_name;
  (void)second_name;
  return true;
}

inline TH1F restore_binning_with_flow(TH1F const& source, TH1F const& reference,
                                      std::string const& context) {
  if (source.GetNbinsX() != reference.GetNbinsX()) {
    throw std::invalid_argument("ChronoSpectra: bin-count mismatch while restoring " + context);
  }
  double underflow = source.GetBinContent(0);
  double underflow_error = source.GetBinError(0);
  double overflow = source.GetBinContent(source.GetNbinsX() + 1);
  double overflow_error = source.GetBinError(source.GetNbinsX() + 1);
  TH1F result = ch::RestoreBinning(source, reference);
  result.SetBinError(0, underflow_error);
  result.SetBinError(result.GetNbinsX() + 1, overflow_error);
  // AddBinContent writes directly to the TH1F storage, including the flow
  // cells.  This also works for a reference histogram whose regular bins have
  // been reset by RestoreBinning.
  result.AddBinContent(0, underflow);
  result.AddBinContent(result.GetNbinsX() + 1, overflow);
  return result;
}

struct GlobPattern {
  std::string bin;
  std::string process;
  std::string parameter;
};

inline std::string normalize_glob_component(std::string component) {
  component = trim(component);
  std::string result;
  for (std::size_t i = 0; i < component.size(); ++i) {
    if (component[i] == '.' && i + 1 < component.size() && component[i + 1] == '*') {
      result.push_back('*');
      ++i;
    } else {
      result.push_back(component[i]);
    }
  }
  return result;
}

inline GlobPattern parse_glob_pattern(std::string pattern) {
  pattern = trim(pattern);
  if (pattern == "all") pattern = "*/*/*";
  std::vector<std::string> pieces = split_trimmed(pattern, '/');
  if (pieces.size() != 3 || pieces[0].empty() || pieces[1].empty() || pieces[2].empty()) {
    throw std::invalid_argument("ChronoSpectra: systematic pattern must have exactly three components: " +
                                quote(pattern));
  }
  GlobPattern result;
  result.bin = normalize_glob_component(pieces[0]);
  result.process = normalize_glob_component(pieces[1]);
  result.parameter = normalize_glob_component(pieces[2]);
  return result;
}

inline bool glob_match(std::string const& pattern, std::string const& value) {
  std::size_t p = 0;
  std::size_t v = 0;
  std::size_t star = std::string::npos;
  std::size_t retry = 0;
  while (v < value.size()) {
    if (p < pattern.size() && pattern[p] == '\\' && p + 1 < pattern.size()) {
      if (pattern[p + 1] == value[v]) { p += 2; ++v; continue; }
    } else if (p < pattern.size() && (pattern[p] == '?' || pattern[p] == value[v])) {
      ++p; ++v; continue;
    } else if (p < pattern.size() && pattern[p] == '*') {
      star = p++;
      retry = v;
      continue;
    }
    if (star == std::string::npos) return false;
    p = star + 1;
    v = ++retry;
  }
  while (p < pattern.size() && pattern[p] == '*') ++p;
  return p == pattern.size();
}

}  // namespace detail
}  // namespace chronospectra
}  // namespace ch

#endif
