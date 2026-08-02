#include "../bin/ChronoSpectraInternal.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace detail = ch::chronospectra::detail;

namespace {

void check(bool condition, std::string const& message) {
  if (!condition) throw std::runtime_error(message);
}

void test_groups() {
  std::set<std::string> names;
  names.insert("binA");
  names.insert("binB");
  names.insert("bin.with.dot");
  detail::NamedGroups groups = detail::parse_named_groups(
      "  region : binA, binA ; dots: .*\\.with\\.dot  ", names, "bin");
  check(groups.groups.at("region").resolved.size() == 1, "group resolution failed");
  check(groups.groups.at("dots").resolved.at(0) == "bin.with.dot", "full-match regex failed");
  check(groups.warnings.size() == 1, "duplicate member warning missing");
  bool failed = false;
  try { detail::parse_named_groups("bad:.*(", names, "bin"); } catch (std::exception const&) { failed = true; }
  check(failed, "invalid regex accepted");
  failed = false;
  try { detail::parse_named_groups("bad:missing", names, "bin"); } catch (std::exception const&) { failed = true; }
  check(failed, "unmatched regex accepted");
  failed = false;
  try { detail::parse_named_groups("bad:/binA", names, "bin"); } catch (std::exception const&) { failed = true; }
  check(failed, "unsafe group name accepted");
}

void test_manifest_and_glob() {
  detail::OutputManifest manifest;
  manifest.add("prefit/binA/total");
  check(manifest.contains("prefit/binA/total"), "manifest did not retain path");
  bool failed = false;
  try { manifest.add("prefit/binA/total"); } catch (std::exception const&) { failed = true; }
  check(failed, "manifest collision accepted");
  failed = false;
  try { manifest.add("prefit/bad name/total"); } catch (std::exception const&) { failed = true; }
  check(failed, "unsafe output component accepted");
  detail::GlobPattern all = detail::parse_glob_pattern("all");
  check(detail::glob_match(all.bin, "binA") && detail::glob_match(all.parameter, "nu"),
        "all glob failed");
  detail::GlobPattern exact = detail::parse_glob_pattern("binA/proc\\*/nu");
  check(detail::glob_match(exact.process, "proc*") && !detail::glob_match(exact.process, "proc1"),
        "escaped glob failed");
  bool malformed = false;
  try { detail::parse_glob_pattern("bin/proc"); } catch (std::exception const&) { malformed = true; }
  check(malformed, "malformed glob accepted");
}

void test_online_correlation() {
  detail::OnlineCorrelation correlation(2);
  correlation.add(std::vector<double>{1., 2.});
  correlation.add(std::vector<double>{2., 4.});
  correlation.add(std::vector<double>{3., 6.});
  std::vector<double> values = correlation.correlation();
  check(std::fabs(values[0] - 1.) < 1e-12, "positive correlation diagonal failed");
  check(std::fabs(values[1] - 1.) < 1e-12, "positive correlation failed");
  check(std::fabs(values[2] - 1.) < 1e-12, "positive correlation failed");
  check(std::fabs(values[3] - 1.) < 1e-12, "positive correlation diagonal failed");

  detail::OnlineCorrelation cancellation(2);
  for (int i = 0; i < 100; ++i) {
    double x = 1.e8 + static_cast<double>(i);
    cancellation.add(std::vector<double>{x, 2.e8 - static_cast<double>(i)});
  }
  values = cancellation.correlation();
  check(values[1] < -0.999999 && values[2] < -0.999999, "cancellation-prone correlation failed");

  detail::OnlineCorrelation zero(2);
  zero.add(std::vector<double>{1., 0.});
  zero.add(std::vector<double>{2., 0.});
  values = zero.correlation();
  check(values[0] == 1. && values[3] == 0. && values[1] == 0. && values[2] == 0.,
        "zero variance policy failed");
  check(detail::OnlineCorrelation(4).correlation().size() == 16, "empty correlation size failed");
}

void test_binning_and_global_state() {
  double edges[] = {0., 0.5, 2., 5.};
  TH1F reference("reference", "reference", 3, edges);
  reference.GetXaxis()->SetBinLabel(1, "low");
  reference.GetXaxis()->SetBinLabel(2, "middle");
  reference.GetXaxis()->SetBinLabel(3, "high");
  TH1F source("source", "source", 3, edges);
  source.SetBinContent(0, 7.);
  source.SetBinError(0, 0.2);
  source.SetBinContent(1, 1.);
  source.SetBinError(1, 0.1);
  source.SetBinContent(2, 2.);
  source.SetBinError(2, 0.2);
  source.SetBinContent(3, 3.);
  source.SetBinError(3, 0.3);
  source.SetBinContent(4, 9.);
  source.SetBinError(4, 0.4);
  source.GetXaxis()->SetBinLabel(1, "low");
  source.GetXaxis()->SetBinLabel(2, "middle");
  source.GetXaxis()->SetBinLabel(3, "high");
  check(source.GetBinContent(0) == 7., "source underflow setup failed: " +
        std::to_string(source.GetBinContent(0)));
  std::string reason;
  check(detail::compatible_binning(reference, source, "reference", "source", &reason),
        "compatible nonuniform binning rejected");
  TH1F restored = detail::restore_binning_with_flow(source, reference, "test");
  check(restored.GetBinContent(0) == 7., "underflow content lost: " +
        std::to_string(restored.GetBinContent(0)));
  check(restored.GetBinContent(4) == 9., "overflow content lost: " +
        std::to_string(restored.GetBinContent(4)));
  check(restored.GetBinError(0) == 0.2 && restored.GetBinError(4) == 0.4, "flow error lost");
  check(std::string(restored.GetXaxis()->GetBinLabel(2)) == "middle", "bin label lost");
  TH1F incompatible("incompatible", "incompatible", 3, 0., 3.);
  check(!detail::compatible_binning(reference, incompatible, "reference", "incompatible", &reason),
        "incompatible edge accepted");
  check(!reason.empty(), "binning failure had no context");

  bool before = TH1::AddDirectoryStatus();
  {
    detail::ScopedAddDirectory guard(!before);
    check(TH1::AddDirectoryStatus() == !before, "AddDirectory guard did not set state");
  }
  check(TH1::AddDirectoryStatus() == before, "AddDirectory guard did not restore state");
}

void assert_parameter(ch::CombineHarvester& cmb, double value, double err_d,
                      double err_u, double range_d, double range_u, bool frozen,
                      std::set<std::string> const& groups, double linked) {
  ch::Parameter const* parameter = cmb.GetParameter("nu");
  check(parameter != nullptr, "parameter disappeared");
  check(std::fabs(parameter->val() - value) < 1e-12 &&
        std::fabs(parameter->err_d() - err_d) < 1e-12 &&
        std::fabs(parameter->err_u() - err_u) < 1e-12,
        "parameter value/error not restored: value=" + std::to_string(parameter->val()) +
        " err_d=" + std::to_string(parameter->err_d()) +
        " err_u=" + std::to_string(parameter->err_u()));
  ch::Parameter *mutable_parameter = const_cast<ch::Parameter*>(parameter);
  check(parameter->range_d() == range_d && parameter->range_u() == range_u &&
        parameter->frozen() == frozen && mutable_parameter->groups() == groups,
        "parameter range/frozen/groups not restored");
  check(mutable_parameter->vars().at(0)->getVal() == linked, "linked RooRealVar not restored");
}

void test_parameter_guard() {
  ch::CombineHarvester cmb;
  cmb.CreateParameterIfEmpty("nu");
  ch::Parameter *parameter = cmb.GetParameter("nu");
  RooRealVar linked("nu", "nu", 1.25);
  parameter->vars().push_back(&linked);
  parameter->set_val(1.25);
  parameter->set_err_d(-0.4);
  parameter->set_err_u(0.7);
  parameter->set_range(-3., 4.);
  parameter->groups().insert("original");
  parameter->set_frozen(false);
  {
    detail::ParameterStateGuard guard(cmb);
    parameter->set_val(9.);
    parameter->set_err_d(-9.);
    parameter->set_err_u(9.);
    parameter->set_range(-9., 9.);
    parameter->groups().clear();
    parameter->groups().insert("changed");
    parameter->set_frozen(true);
    guard.restore();
  }
  std::set<std::string> original_groups;
  original_groups.insert("original");
  assert_parameter(cmb, 1.25, -0.4, 0.7, -3., 4., false,
                   original_groups, 1.25);

  parameter->set_frozen(false);
  parameter->set_val(2.5);
  linked.setVal(2.5);
  parameter->set_frozen(true);
  std::set<std::string> frozen_groups;
  frozen_groups.insert("frozen");
  parameter->groups() = frozen_groups;
  bool threw = false;
  try {
    detail::ParameterStateGuard guard(cmb);
    parameter->set_frozen(false);
    parameter->set_val(-4.);
    parameter->groups().clear();
    throw std::runtime_error("deliberate test exception");
  } catch (std::exception const&) {
    threw = true;
  }
  check(threw, "deliberate exception did not execute");
  assert_parameter(cmb, 2.5, -0.4, 0.7, -3., 4., true, frozen_groups, 2.5);
}

}  // namespace

int main() {
  try {
    test_groups();
    test_manifest_and_glob();
    test_online_correlation();
    test_binning_and_global_state();
    test_parameter_guard();
    std::cout << "ChronoSpectraInternalTest: all invariant tests passed\n";
    return 0;
  } catch (std::exception const& error) {
    std::cerr << "ChronoSpectraInternalTest: FAIL: " << error.what() << "\n";
    return 1;
  }
}
