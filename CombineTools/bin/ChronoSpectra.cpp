// ChronoSpectra concept and original implementation: Mohammad Abrar Wadud (2024).
//
// This executable is a narrow composition layer over the existing
// CombineHarvester evaluation API.  In particular, it intentionally does not
// duplicate the target's RooFit, uncertainty, or two-dimensional algorithms.

#include "ChronoSpectraInternal.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooSimultaneous.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMatrixDSym.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "CombineHarvester/CombineTools/interface/ParseCombineWorkspace.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace detail = ch::chronospectra::detail;

namespace {

struct Config {
  std::string workspace;
  std::string datacard;
  std::string output;
  std::string dataset;
  std::string mass;
  std::string fitresult;
  unsigned samples;
  bool postfit;
  bool skipprefit;
  std::string freeze;
  std::string group_bins;
  std::string group_procs;
  bool skip_obs;
  bool get_rate_corr;
  bool get_hist_bin_corr;
  bool sep_proc_hists;
  bool sep_bin_hists;
  bool sep_proc_hist_bin_corr;
  bool sep_bin_hist_bin_corr;
  bool sep_bin_rate_corr;
  bool store_syst;
  std::string plot_syst;
  std::string syst_save_dir;
  bool logy;
  std::string log_level;

  Config()
      : dataset("data_obs"),
        mass("125"),
        samples(2000),
        postfit(false),
        skipprefit(false),
        skip_obs(false),
        get_rate_corr(true),
        get_hist_bin_corr(true),
        sep_proc_hists(false),
        sep_bin_hists(false),
        sep_proc_hist_bin_corr(false),
        sep_bin_hist_bin_corr(false),
        sep_bin_rate_corr(false),
        store_syst(false),
        syst_save_dir("shapeSystPlots"),
        logy(false),
        log_level("info") {}
};

struct Selection {
  std::string name;
  std::vector<std::string> bins;
  bool is_bin_group;
  bool is_actual_bin;
};

struct HistogramProduct {
  std::string path;
  TH1F histogram;
};

struct MatrixProduct {
  std::string path;
  TH2F matrix;
};

struct SystematicProduct {
  std::string bin;
  std::string process;
  std::string parameter;
  TH1F nominal;
  TH1F up;
  TH1F down;
};

struct PreparedModel {
  std::unique_ptr<TFile> workspace_file;
  std::unique_ptr<RooFitResult> fit;
  ch::CombineHarvester workspace;
  ch::CombineHarvester card;
  std::map<std::string, TH1F> references;
  detail::NamedGroups bin_groups;
  detail::NamedGroups process_groups;
  std::vector<Selection> selections;
  std::set<std::string> grouped_bins;
  std::set<std::string> grouped_processes;
};

struct ParsedArguments {
  Config config;
  po::options_description description;
  po::variables_map variables;
  std::string samples_text;
  bool help;

  ParsedArguments()
      : description("ChronoSpectra options"), samples_text("2000"), help(false) {}
};

void add_bool_option(po::options_description& options, char const* name, bool* target,
                     char const* description, bool default_value) {
  options.add_options()(name, po::value<bool>(target)->default_value(default_value)
                                  ->implicit_value(true), description);
}

ParsedArguments parse_arguments(int argc, char** argv) {
  ParsedArguments result;
  Config& cfg = result.config;
  result.description.add_options()
      ("help,h", "Print this help and exit")
      ("workspace,w", po::value<std::string>(&cfg.workspace)->required(),
       "ROOT workspace file containing object w")
      ("datacard,d", po::value<std::string>(&cfg.datacard)->required(),
       "original datacard used for physical binning")
      ("output,o", po::value<std::string>(&cfg.output)->required(),
       "ROOT output file, committed only after success")
      ("dataset", po::value<std::string>(&cfg.dataset)->default_value(cfg.dataset),
       "workspace dataset and output key (default: data_obs)")
      ("mass,m", po::value<std::string>(&cfg.mass)->default_value(cfg.mass),
       "datacard mass (default: 125)")
      ("fitresult,f", po::value<std::string>(&cfg.fitresult)->default_value(std::string()),
       "fit ROOT file and object, file.root:object")
      ("samples", po::value<std::string>(&result.samples_text)
                    ->default_value(std::to_string(cfg.samples)),
       "post-fit samples; zero selects target unsampled errors (default: 2000)")
      ("seed", po::value<std::string>(), "RooFit random seed")
      ("freeze", po::value<std::string>(&cfg.freeze)->default_value(std::string()),
       "comma-separated PARAM or PARAM=VALUE freezes")
      ("groupBins", po::value<std::string>(&cfg.group_bins)->default_value(std::string()),
       "named bin groups: name:member,member;name:member")
      ("groupProcs", po::value<std::string>(&cfg.group_procs)->default_value(std::string()),
       "named process groups: name:member,member;name:member")
      ("plotSyst", po::value<std::string>(&cfg.plot_syst)->default_value(std::string()),
       "systematic plot patterns (requires --storeSyst)")
      ("systSaveDir", po::value<std::string>(&cfg.syst_save_dir)->default_value(cfg.syst_save_dir),
       "systematic PNG destination (default: shapeSystPlots)");
  add_bool_option(result.description, "postfit", &cfg.postfit,
                  "also emit post-fit products", false);
  add_bool_option(result.description, "skipprefit", &cfg.skipprefit,
                  "omit pre-fit products; requires --postfit", false);
  add_bool_option(result.description, "skipObs", &cfg.skip_obs,
                  "replace observed data with total pseudo-data", false);
  add_bool_option(result.description, "getRateCorr", &cfg.get_rate_corr,
                  "emit target rate correlations when sampled", true);
  add_bool_option(result.description, "getHistBinCorr", &cfg.get_hist_bin_corr,
                  "emit local histogram-bin correlations when sampled", true);
  add_bool_option(result.description, "sepProcHists", &cfg.sep_proc_hists,
                  "restore individual process histograms suppressed by groups", false);
  add_bool_option(result.description, "sepBinHists", &cfg.sep_bin_hists,
                  "restore actual-bin histograms suppressed by groups", false);
  add_bool_option(result.description, "sepProcHistBinCorr", &cfg.sep_proc_hist_bin_corr,
                  "restore grouped individual-process bin matrices", false);
  add_bool_option(result.description, "sepBinHistBinCorr", &cfg.sep_bin_hist_bin_corr,
                  "restore grouped actual-bin bin matrices", false);
  add_bool_option(result.description, "sepBinRateCorr", &cfg.sep_bin_rate_corr,
                  "restore grouped actual-bin rate matrices", false);
  add_bool_option(result.description, "storeSyst", &cfg.store_syst,
                  "store pre-fit systematic variations", false);
  add_bool_option(result.description, "logy", &cfg.logy,
                  "use logarithmic systematic-plot scale", false);
  result.description.add_options()
      ("logLevel", po::value<std::string>(&cfg.log_level)->default_value(cfg.log_level),
       "log level: info, warn, or error")
      ("log-level", po::value<std::string>(&cfg.log_level),
       "compatibility spelling of --logLevel");

  po::variables_map variables;
  po::parsed_options parsed = po::command_line_parser(argc, argv)
                                  .options(result.description).run();
  po::store(parsed, variables);
  if (variables.count("help")) {
    result.help = true;
    return result;
  }
  po::notify(variables);
  cfg.samples = detail::parse_unsigned(result.samples_text, "samples");
  if (variables.count("seed")) {
    detail::parse_unsigned(variables.at("seed").as<std::string>(), "seed");
  }
  result.variables = variables;
  return result;
}

void validate_config(Config const& cfg) {
  if (cfg.workspace.empty() || cfg.datacard.empty() || cfg.output.empty()) {
    throw std::invalid_argument("ChronoSpectra: workspace, datacard, and output are required");
  }
  if (cfg.skipprefit && !cfg.postfit) {
    throw std::invalid_argument("ChronoSpectra: --skipprefit requires --postfit");
  }
  if (cfg.postfit && cfg.fitresult.empty()) {
    throw std::invalid_argument("ChronoSpectra: --postfit requires --fitresult=file.root:object");
  }
  if (!cfg.plot_syst.empty() && !cfg.store_syst) {
    throw std::invalid_argument("ChronoSpectra: --plotSyst requires --storeSyst");
  }
  if (cfg.store_syst && cfg.skipprefit) {
    throw std::invalid_argument("ChronoSpectra: --storeSyst requires pre-fit products");
  }
  if (cfg.log_level != "info" && cfg.log_level != "warn" && cfg.log_level != "error") {
    throw std::invalid_argument("ChronoSpectra: unknown log level " + detail::quote(cfg.log_level));
  }
  detail::require_safe_path_component(cfg.dataset, "dataset");
  if (cfg.dataset == "signal" || cfg.dataset == "background" || cfg.dataset == "total") {
    throw std::invalid_argument("ChronoSpectra: dataset collides with an aggregate key");
  }
  if (!cfg.plot_syst.empty()) {
    std::vector<std::string> patterns = detail::split_trimmed(cfg.plot_syst, ',');
    for (std::vector<std::string>::const_iterator it = patterns.begin(); it != patterns.end(); ++it) {
      detail::parse_glob_pattern(*it);
    }
  }
}

std::set<std::string> set_difference_names(std::set<std::string> const& left,
                                           std::set<std::string> const& right) {
  std::set<std::string> result;
  std::set_difference(left.begin(), left.end(), right.begin(), right.end(),
                      std::inserter(result, result.end()));
  return result;
}

void validate_name_collisions(Config const& cfg, PreparedModel& model) {
  std::set<std::string> reserved;
  reserved.insert("signal");
  reserved.insert("background");
  reserved.insert("total");
  reserved.insert(cfg.dataset);
  std::set<std::string> all_names = model.workspace.bin_set();
  std::set<std::string> processes = model.workspace.process_set();
  all_names.insert(processes.begin(), processes.end());
  for (std::map<std::string, detail::NamedGroup>::const_iterator it = model.bin_groups.groups.begin();
       it != model.bin_groups.groups.end(); ++it) all_names.insert(it->first);
  for (std::map<std::string, detail::NamedGroup>::const_iterator it = model.process_groups.groups.begin();
       it != model.process_groups.groups.end(); ++it) all_names.insert(it->first);
  for (std::set<std::string>::const_iterator it = all_names.begin(); it != all_names.end(); ++it) {
    detail::require_safe_path_component(*it, "model");
  }
  std::set<std::string> group_names;
  for (std::map<std::string, detail::NamedGroup>::const_iterator it = model.bin_groups.groups.begin();
       it != model.bin_groups.groups.end(); ++it) group_names.insert(it->first);
  for (std::map<std::string, detail::NamedGroup>::const_iterator it = model.process_groups.groups.begin();
       it != model.process_groups.groups.end(); ++it) {
    if (!group_names.insert(it->first).second) {
      throw std::invalid_argument("ChronoSpectra: bin/process group name collision " + detail::quote(it->first));
    }
  }
  for (std::set<std::string>::const_iterator it = group_names.begin(); it != group_names.end(); ++it) {
    if (reserved.count(*it)) {
      throw std::invalid_argument("ChronoSpectra: group name collides with reserved output key " +
                                  detail::quote(*it));
    }
    if (model.workspace.bin_set().count(*it) || processes.count(*it)) {
      throw std::invalid_argument("ChronoSpectra: group name collides with an actual model name " +
                                  detail::quote(*it));
    }
  }
}

std::vector<std::string> intersection(std::vector<std::string> const& available,
                                      std::vector<std::string> const& requested) {
  std::set<std::string> requested_set(requested.begin(), requested.end());
  std::vector<std::string> result;
  for (std::vector<std::string>::const_iterator it = available.begin(); it != available.end(); ++it) {
    if (requested_set.count(*it)) result.push_back(*it);
  }
  return result;
}

std::vector<std::string> as_vector(std::set<std::string> const& names) {
  return std::vector<std::string>(names.begin(), names.end());
}

ch::CombineHarvester select(ch::CombineHarvester& cmb,
                            std::vector<std::string> const& bins,
                            std::vector<std::string> const& processes) {
  return cmb.cp().bin(bins).process(processes);
}

TH1F restore_shape(TH1F shape, TH1F const& reference, std::string const& context) {
  for (int i = 1; i <= shape.GetNbinsX(); ++i) {
    if (!std::isfinite(shape.GetBinContent(i)) || !std::isfinite(shape.GetBinError(i))) {
      throw std::runtime_error("ChronoSpectra: non-finite shape content/error in " + context);
    }
  }
  return detail::restore_binning_with_flow(shape, reference, context);
}

TH1F pseudo_observation(TH1F const& total, std::string const& dataset,
                        std::string const& context) {
  TH1F result = total;
  result.SetTitle(dataset.c_str());
  for (int i = 1; i <= result.GetNbinsX(); ++i) {
    if (!std::isfinite(result.GetBinContent(i)) || result.GetBinContent(i) < 0.) {
      throw std::runtime_error("ChronoSpectra: cannot construct Poisson pseudo-data with negative/non-finite " +
                               context + " bin " + std::to_string(i));
    }
    result.SetBinError(i, 0.);
  }
  double total_rate = result.Integral();
  if (!std::isfinite(total_rate) || total_rate < 0.) {
    throw std::runtime_error("ChronoSpectra: negative/non-finite pseudo-data total in " + context);
  }
  result.SetBinErrorOption(TH1::kPoisson);
  result.SetBinContent(0, std::sqrt(total_rate));
  result.SetBinError(0, 0.);
  return result;
}

void add_histogram(detail::OutputManifest& manifest,
                   std::vector<HistogramProduct>& products,
                   std::string const& path, TH1F histogram, std::string const& title) {
  manifest.add(path);
  histogram.SetTitle(title.c_str());
  products.push_back(HistogramProduct{path, histogram});
}

void add_matrix(detail::OutputManifest& manifest,
                std::vector<MatrixProduct>& products,
                std::string const& path, TH2F matrix, std::string const& title) {
  manifest.add(path);
  matrix.SetTitle(title.c_str());
  products.push_back(MatrixProduct{path, matrix});
}

bool meaningfully_different(TH1F const& first, TH1F const& second) {
  if (first.GetNbinsX() != second.GetNbinsX()) return true;
  for (int i = 1; i <= first.GetNbinsX(); ++i) {
    double scale = std::max(1., std::max(std::fabs(first.GetBinContent(i)),
                                         std::fabs(second.GetBinContent(i))));
    if (std::fabs(first.GetBinContent(i) - second.GetBinContent(i)) > 1e-12 * scale) return true;
  }
  return false;
}

void store_systematics(Config const& cfg, PreparedModel& model,
                       detail::OutputManifest& manifest,
                       std::vector<HistogramProduct>& histograms,
                       std::vector<SystematicProduct>& systematics) {
  (void)cfg;
  std::vector<ch::Parameter> parameters = model.workspace.GetParameters();
  for (std::map<std::string, TH1F>::const_iterator reference = model.references.begin();
       reference != model.references.end(); ++reference) {
    std::vector<std::string> one_bin(1, reference->first);
    ch::CombineHarvester bin_scope = select(model.workspace, one_bin,
                                            as_vector(model.workspace.process_set()));
    std::set<std::string> bin_processes = bin_scope.process_set();
    for (std::set<std::string>::const_iterator process = bin_processes.begin();
         process != bin_processes.end(); ++process) {
      std::vector<std::string> one_process(1, *process);
      ch::CombineHarvester process_scope = select(model.workspace, one_bin, one_process);
      std::string base = "systematics/" + reference->first + "/" + *process;
      TH1F nominal = restore_shape(process_scope.GetShapeWithUncertainty(), reference->second,
                                   "systematic nominal " + reference->first + "/" + *process);
      double nominal_uncertainty = process_scope.GetUncertainty();
      nominal.AddBinContent(0, nominal_uncertainty - nominal.GetBinContent(0));
      nominal.SetBinError(0, 0.);
      nominal.SetTitle(process->c_str());
      add_histogram(manifest, histograms, base, nominal, *process);

      for (std::vector<ch::Parameter>::const_iterator saved = parameters.begin();
           saved != parameters.end(); ++saved) {
        ch::Parameter* parameter = model.workspace.GetParameter(saved->name());
        if (!parameter || parameter->frozen()) continue;
        if (!std::isfinite(parameter->val()) || !std::isfinite(parameter->err_d()) ||
            !std::isfinite(parameter->err_u())) {
          throw std::runtime_error("ChronoSpectra: non-finite systematic state for parameter " + saved->name());
        }
        parameter->set_val(saved->val() + saved->err_u());
        TH1F up = restore_shape(process_scope.GetShape(), reference->second,
                                "systematic up " + reference->first + "/" + *process + "/" + saved->name());
        parameter->set_val(saved->val() + saved->err_d());
        TH1F down = restore_shape(process_scope.GetShape(), reference->second,
                                  "systematic down " + reference->first + "/" + *process + "/" + saved->name());
        parameter->set_val(saved->val());
        if (!meaningfully_different(nominal, up) && !meaningfully_different(nominal, down)) continue;
        detail::require_safe_path_component(saved->name(), "parameter");
        std::string path_base = base + "_syst/" + saved->name();
        up.SetTitle((saved->name() + " Up").c_str());
        down.SetTitle((saved->name() + " Down").c_str());
        add_histogram(manifest, histograms, path_base + "_Up", up, saved->name() + " Up");
        add_histogram(manifest, histograms, path_base + "_Down", down, saved->name() + " Down");
        systematics.push_back(SystematicProduct{reference->first, *process, saved->name(), nominal, up, down});
      }
    }
  }
}

void plot_systematics(Config const& cfg, std::vector<SystematicProduct> const& systematics) {
  if (cfg.plot_syst.empty()) return;
  std::vector<detail::GlobPattern> patterns;
  std::vector<std::string> pattern_strings = detail::split_trimmed(cfg.plot_syst, ',');
  for (std::vector<std::string>::const_iterator pattern = pattern_strings.begin();
       pattern != pattern_strings.end(); ++pattern) patterns.push_back(detail::parse_glob_pattern(*pattern));

  std::vector<SystematicProduct const*> selected;
  std::set<std::string> matched_patterns;
  for (std::vector<SystematicProduct>::const_iterator systematic = systematics.begin();
       systematic != systematics.end(); ++systematic) {
    bool matches = false;
    for (std::size_t i = 0; i < patterns.size(); ++i) {
      if (detail::glob_match(patterns[i].bin, systematic->bin) &&
          detail::glob_match(patterns[i].process, systematic->process) &&
          detail::glob_match(patterns[i].parameter, systematic->parameter)) {
        matches = true;
        matched_patterns.insert(pattern_strings[i]);
      }
    }
    if (matches) selected.push_back(&*systematic);
  }
  for (std::vector<std::string>::const_iterator pattern = pattern_strings.begin();
       pattern != pattern_strings.end(); ++pattern) {
    if (!matched_patterns.count(*pattern)) {
      std::cerr << "ChronoSpectra: warning: systematic plot pattern matched no stored variation: "
                << detail::quote(*pattern) << "\n";
    }
  }
  if (selected.empty()) return;

  fs::path destination(cfg.syst_save_dir);
  if (destination.empty()) throw std::runtime_error("ChronoSpectra: empty systematic plot directory");
  if (fs::exists(destination) && !fs::is_directory(destination)) {
    throw std::runtime_error("ChronoSpectra: systematic plot destination is not a directory: " + destination.string());
  }
  fs::create_directories(destination);
  fs::path staging = destination.parent_path() / fs::unique_path(destination.filename().string() + ".chrono-%%%%%%.staging");
  fs::create_directories(staging);
  std::vector<fs::path> final_paths;
  std::vector<fs::path> moved_paths;
  try {
    for (std::vector<SystematicProduct const*>::const_iterator systematic = selected.begin();
         systematic != selected.end(); ++systematic) {
      std::string basename = (*systematic)->bin + "__" + (*systematic)->process + "__" +
                             (*systematic)->parameter + ".png";
      fs::path final = destination / basename;
      if (fs::exists(final)) throw std::runtime_error("ChronoSpectra: systematic plot already exists: " + final.string());
      final_paths.push_back(final);
      fs::path staged = staging / basename;
      TCanvas canvas("ChronoSpectraSystematic", "ChronoSpectra systematic", 2800, 2400);
      TPad upper("upper", "upper", 0., 0.40, 1., 1.);
      TPad lower("lower", "lower", 0., 0., 1., 0.40);
      upper.SetLeftMargin(0.25);
      upper.SetRightMargin(0.05);
      upper.SetBottomMargin(0.015);
      upper.SetTopMargin(0.10);
      upper.SetGrid(1, 1);
      lower.SetLeftMargin(0.25);
      lower.SetRightMargin(0.05);
      lower.SetBottomMargin(0.38);
      lower.SetTopMargin(0.0);
      lower.SetGrid(1, 1);
      upper.Draw();
      lower.Draw();
      upper.cd();
      TH1F nominal = (*systematic)->nominal;
      TH1F up = (*systematic)->up;
      TH1F down = (*systematic)->down;
      nominal.SetStats(false);
      up.SetStats(false);
      down.SetStats(false);
      std::string plot_title = (*systematic)->bin + " / " + (*systematic)->process +
                               " / " + (*systematic)->parameter;
      nominal.SetTitle(plot_title.c_str());
      nominal.GetXaxis()->SetLabelSize(0.);
      nominal.GetXaxis()->SetTitleSize(0.);
      nominal.GetYaxis()->SetTitle("Events");
      nominal.GetYaxis()->CenterTitle();
      nominal.GetYaxis()->SetTitleOffset(0.96);
      nominal.GetYaxis()->SetTitleSize(0.10);
      nominal.GetYaxis()->SetLabelSize(0.085);
      nominal.GetYaxis()->SetMoreLogLabels();
      nominal.SetLineColor(kBlack);
      nominal.SetLineWidth(5);
      up.SetLineColor(kRed + 1);
      up.SetLineWidth(5);
      down.SetLineColor(kBlue + 1);
      down.SetLineWidth(5);

      double y_min = std::numeric_limits<double>::max();
      double y_max = -std::numeric_limits<double>::max();
      double positive_min = std::numeric_limits<double>::max();
      bool all_positive = true;
      std::vector<TH1F const*> shapes;
      shapes.push_back(&nominal);
      shapes.push_back(&up);
      shapes.push_back(&down);
      for (std::vector<TH1F const*>::const_iterator shape = shapes.begin();
           shape != shapes.end(); ++shape) {
        for (int bin = 1; bin <= (*shape)->GetNbinsX(); ++bin) {
          double value = (*shape)->GetBinContent(bin);
          if (!std::isfinite(value)) continue;
          y_min = std::min(y_min, value);
          y_max = std::max(y_max, value);
          if (value > 0.) {
            positive_min = std::min(positive_min, value);
          } else {
            all_positive = false;
          }
        }
      }
      if (y_min == std::numeric_limits<double>::max() ||
          y_max == -std::numeric_limits<double>::max()) {
        throw std::runtime_error("ChronoSpectra: systematic plot has no finite regular-bin values");
      }

      bool use_logy = cfg.logy && all_positive &&
                      positive_min != std::numeric_limits<double>::max();
      if (cfg.logy && !use_logy) {
        std::cerr << "ChronoSpectra: warning: non-positive systematic-plot content in "
                  << plot_title << "; using linear y scale\n";
      }
      if (use_logy) {
        upper.SetLogy(true);
        nominal.SetMinimum(std::max(0.8 * positive_min, positive_min / 10.));
        nominal.SetMaximum(std::max(1.2 * y_max, 2. * positive_min));
      } else {
        double span = y_max - y_min;
        double padding = span > 0. ? 0.05 * span : std::max(1., std::fabs(y_max) * 0.10);
        nominal.SetMinimum(y_min - padding);
        nominal.SetMaximum(y_max + padding);
      }

      nominal.Draw("hist");
      up.Draw("hist same");
      down.Draw("hist same");

      auto format_value = [](double value) {
        if (!std::isfinite(value)) return std::string("n/a");
        std::ostringstream formatted;
        double magnitude = std::fabs(value);
        if ((magnitude > 0. && magnitude < 0.01) || magnitude >= 1000.) {
          formatted << std::scientific << std::setprecision(2) << value;
        } else {
          formatted << std::fixed << std::setprecision(magnitude >= 10. ? 1 : 2) << value;
        }
        return formatted.str();
      };
      auto format_change = [&format_value](double value, double nominal_value) {
        if (!std::isfinite(value) || !std::isfinite(nominal_value) || nominal_value == 0.) {
          return std::string("n/a");
        }
        double change = 100. * (value - nominal_value) / nominal_value;
        return std::string(change >= 0. ? "+" : "") + format_value(change) + "%";
      };
      double nominal_integral = nominal.Integral();
      double up_integral = up.Integral();
      double down_integral = down.Integral();
      TLegend legend(0.57, 0.67, 0.95, 0.90, "", "NBNDC");
      legend.SetBorderSize(0);
      legend.SetFillStyle(0);
      legend.SetTextSize(0.046);
      legend.AddEntry(&nominal, ("Nominal (n=" + format_value(nominal_integral) + ")").c_str(), "l");
      legend.AddEntry(&up, ("Up (n=" + format_value(up_integral) + ", " +
                            format_change(up_integral, nominal_integral) + ")").c_str(), "l");
      legend.AddEntry(&down, ("Down (n=" + format_value(down_integral) + ", " +
                              format_change(down_integral, nominal_integral) + ")").c_str(), "l");
      legend.Draw();

      lower.cd();
      TH1F relative_up = up;
      TH1F relative_down = down;
      relative_up.Reset();
      relative_down.Reset();
      relative_up.SetStats(false);
      relative_down.SetStats(false);
      relative_up.SetTitle("");
      relative_up.GetYaxis()->SetTitle("Variation (%)");
      relative_up.GetYaxis()->CenterTitle();
      relative_up.GetYaxis()->SetTitleOffset(0.67);
      relative_up.GetYaxis()->SetTitleSize(0.145);
      relative_up.GetYaxis()->SetLabelSize(0.13);
      relative_up.GetYaxis()->SetNdivisions(505);
      relative_up.GetYaxis()->SetMaxDigits(3);
      relative_up.GetXaxis()->CenterTitle();
      relative_up.GetXaxis()->SetTitleOffset(1.0);
      relative_up.GetXaxis()->SetTitleSize(0.145);
      relative_up.GetXaxis()->SetLabelSize(0.13);
      relative_up.SetLineColor(kRed + 1);
      relative_up.SetLineWidth(5);
      relative_down.SetLineColor(kBlue + 1);
      relative_down.SetLineWidth(5);
      bool undefined_relative_bin = false;
      for (int bin = 1; bin <= relative_up.GetNbinsX(); ++bin) {
        double nominal_value = nominal.GetBinContent(bin);
        if (std::isfinite(nominal_value) && nominal_value > 0.) {
          relative_up.SetBinContent(
              bin, 100. * (up.GetBinContent(bin) - nominal_value) / nominal_value);
          relative_down.SetBinContent(
              bin, 100. * (down.GetBinContent(bin) - nominal_value) / nominal_value);
        } else {
          relative_up.SetBinContent(bin, 0.);
          relative_down.SetBinContent(bin, 0.);
          undefined_relative_bin = true;
        }
      }
      if (undefined_relative_bin) {
        std::cerr << "ChronoSpectra: warning: non-positive nominal bin in relative systematic plot "
                  << plot_title << "; that bin is shown at zero\n";
      }
      double relative_min = std::numeric_limits<double>::max();
      double relative_max = -std::numeric_limits<double>::max();
      for (TH1F const* relative : {&relative_up, &relative_down}) {
        for (int bin = 1; bin <= relative->GetNbinsX(); ++bin) {
          relative_min = std::min(relative_min, relative->GetBinContent(bin));
          relative_max = std::max(relative_max, relative->GetBinContent(bin));
        }
      }
      double relative_extent = std::max(std::fabs(relative_min), std::fabs(relative_max));
      double relative_padding = relative_extent > 0. ? 0.20 * relative_extent : 1.;
      relative_up.SetMinimum(std::min(-relative_padding, relative_min - relative_padding));
      relative_up.SetMaximum(std::max(relative_padding, relative_max + relative_padding));
      double x_min = relative_up.GetXaxis()->GetXmin();
      double x_max = relative_up.GetXaxis()->GetXmax();
      TLine zero_line(x_min, 0., x_max, 0.);
      zero_line.SetLineColor(kBlack);
      zero_line.SetLineWidth(4);
      relative_up.Draw("hist");
      relative_down.Draw("hist same");
      zero_line.Draw("same");
      upper.RedrawAxis();
      lower.RedrawAxis();
      canvas.Print(staged.string().c_str(), "png");
      if (!fs::exists(staged) || fs::file_size(staged) == 0) throw std::runtime_error("ChronoSpectra: failed to render systematic plot " + staged.string());
    }
    for (std::size_t i = 0; i < final_paths.size(); ++i) {
      fs::rename(staging / final_paths[i].filename(), final_paths[i]);
      moved_paths.push_back(final_paths[i]);
    }
    fs::remove_all(staging);
  } catch (...) {
    for (std::vector<fs::path>::const_iterator moved = moved_paths.begin(); moved != moved_paths.end(); ++moved) {
      if (fs::exists(*moved)) fs::remove(*moved);
    }
    if (fs::exists(staging)) fs::remove_all(staging);
    throw;
  }
}

TH2F rate_correlation(ch::CombineHarvester& selection,
                      RooFitResult const& fit, unsigned samples,
                      std::string const& context) {
  TH2F result = selection.GetRateCorrelation(fit, samples);
  bool sanitized = false;
  for (int x = 1; x <= result.GetNbinsX(); ++x) {
    for (int y = 1; y <= result.GetNbinsY(); ++y) {
      if (!std::isfinite(result.GetBinContent(x, y))) {
        result.SetBinContent(x, y, 0.);
        sanitized = true;
      }
    }
  }
  if (sanitized) {
    std::cerr << "ChronoSpectra: warning: non-finite target rate-correlation cells replaced with zero in "
              << context << "\n";
  }
  return result;
}

TH2F parameter_correlation(RooFitResult const& fit) {
  RooArgList const& parameters = fit.floatParsFinal();
  int n = parameters.getSize();
  TH2F result("parCorrMat", "Parameter Correlation Matrix", n, 0.5, n + 0.5,
              n, 0.5, n + 0.5);
  TMatrixDSym correlation = fit.correlationMatrix();
  for (int i = 0; i < n; ++i) {
    RooAbsArg const* parameter = parameters.at(i);
    result.GetXaxis()->SetBinLabel(i + 1, parameter->GetName());
    result.GetYaxis()->SetBinLabel(i + 1, parameter->GetName());
    for (int j = 0; j < n; ++j) result.SetBinContent(i + 1, j + 1, correlation(i, j));
  }
  return result;
}

TH2F histogram_bin_correlation(ch::CombineHarvester& selection,
                               RooFitResult const& fit, unsigned samples,
                               std::string const& context) {
  TH1F nominal = selection.GetShape();
  int bins = nominal.GetNbinsX();
  detail::OnlineCorrelation accumulator(static_cast<std::size_t>(bins));
  {
    detail::ParameterStateGuard sampling_state(selection);
    RooArgList const& random_parameters = fit.randomizePars();
    std::vector<RooRealVar const*> random_variables;
    std::vector<ch::Parameter*> target_parameters;
    for (int i = 0; i < random_parameters.getSize(); ++i) {
      RooRealVar const* variable = dynamic_cast<RooRealVar const*>(random_parameters.at(i));
      if (!variable) continue;
      random_variables.push_back(variable);
      target_parameters.push_back(selection.GetParameter(variable->GetName()));
    }
    for (unsigned sample = 0; sample < samples; ++sample) {
      fit.randomizePars();
      for (std::size_t i = 0; i < random_variables.size(); ++i) {
        if (target_parameters[i]) target_parameters[i]->set_val(random_variables[i]->getVal());
      }
      TH1F sampled = selection.GetShape();
      if (sampled.GetNbinsX() != bins) {
        throw std::runtime_error("ChronoSpectra: sampled shape bin count changed in " + context);
      }
      std::vector<double> values(static_cast<std::size_t>(bins));
      for (int bin = 1; bin <= bins; ++bin) {
        values[static_cast<std::size_t>(bin - 1)] = sampled.GetBinContent(bin);
      }
      accumulator.add(values);
    }
    sampling_state.restore();
  }
  TH2F result("HistBinCorr", "Histogram Bin Correlation Matrix", bins, 0., bins, bins, 0., bins);
  std::vector<double> correlations = accumulator.correlation();
  for (int i = 0; i < bins; ++i) {
    std::string label = "Bin " + std::to_string(i + 1);
    result.GetXaxis()->SetBinLabel(i + 1, label.c_str());
    result.GetYaxis()->SetBinLabel(i + 1, label.c_str());
    for (int j = 0; j < bins; ++j) result.SetBinContent(i + 1, j + 1,
                                                         correlations[static_cast<std::size_t>(i * bins + j)]);
  }
  return result;
}

std::vector<std::string> parse_freeze_tokens(std::string const& definition,
                                             ch::CombineHarvester& cmb) {
  std::vector<std::string> names;
  if (detail::trim(definition).empty()) return names;
  std::set<std::string> seen;
  std::vector<std::string> entries = detail::split_trimmed(definition, ',');
  for (std::vector<std::string>::const_iterator it = entries.begin(); it != entries.end(); ++it) {
    if (it->empty()) throw std::invalid_argument("ChronoSpectra: empty freeze entry");
    std::size_t equals = it->find('=');
    if (equals != it->rfind('=')) throw std::invalid_argument("ChronoSpectra: malformed freeze entry " + detail::quote(*it));
    std::string name = detail::trim(equals == std::string::npos ? *it : it->substr(0, equals));
    if (name.empty()) throw std::invalid_argument("ChronoSpectra: empty freeze parameter name");
    if (!seen.insert(name).second) throw std::invalid_argument("ChronoSpectra: duplicate freeze parameter " + detail::quote(name));
    ch::Parameter* parameter = cmb.GetParameter(name);
    if (!parameter) throw std::invalid_argument("ChronoSpectra: unknown freeze parameter " + detail::quote(name));
    if (equals != std::string::npos) {
      double value = detail::parse_finite_double(it->substr(equals + 1), "freeze " + name);
      parameter->set_frozen(false);
      parameter->set_val(value);
    }
    parameter->set_frozen(true);
    names.push_back(name);
  }
  return names;
}

void validate_fit(RooFitResult const& fit, std::string const& source) {
  if (fit.floatParsFinal().getSize() == 0) {
    throw std::invalid_argument("ChronoSpectra: fit result has no floating final parameters: " + source);
  }
  TMatrixDSym correlation = fit.correlationMatrix();
  int n = fit.floatParsFinal().getSize();
  if (correlation.GetNrows() != n || correlation.GetNcols() != n) {
    throw std::invalid_argument("ChronoSpectra: fit correlation dimensions do not match floating parameters: " + source);
  }
}

PreparedModel prepare_model(Config const& cfg) {
  PreparedModel model;
  model.workspace_file.reset(new TFile(cfg.workspace.c_str(), "READ"));
  if (!model.workspace_file->IsOpen() || model.workspace_file->IsZombie()) {
    throw std::runtime_error("ChronoSpectra: cannot open workspace file " + cfg.workspace);
  }
  RooWorkspace* workspace = dynamic_cast<RooWorkspace*>(model.workspace_file->Get("w"));
  if (!workspace) throw std::runtime_error("ChronoSpectra: workspace object 'w' missing from " + cfg.workspace);
  RooStats::ModelConfig* model_config = dynamic_cast<RooStats::ModelConfig*>(workspace->genobj("ModelConfig"));
  if (!model_config) throw std::runtime_error("ChronoSpectra: ModelConfig missing from " + cfg.workspace);
  if (!dynamic_cast<RooSimultaneous*>(model_config->GetPdf())) {
    throw std::runtime_error("ChronoSpectra: ModelConfig PDF is not a RooSimultaneous");
  }
  if (!workspace->data(cfg.dataset.c_str())) {
    throw std::runtime_error("ChronoSpectra: dataset " + detail::quote(cfg.dataset) + " missing from workspace");
  }

  model.workspace.SetFlag("workspaces-use-clone", true);
  ch::ParseCombineWorkspace(model.workspace, *workspace, "ModelConfig", cfg.dataset, false);
  model.workspace.FilterProcs([](ch::Process* process) {
    return !process->shape() && !process->data() && !process->pdf();
  });
  std::set<std::string> bins = model.workspace.bin_set();
  std::set<std::string> processes = model.workspace.process_set();
  if (bins.empty() || processes.empty()) throw std::runtime_error("ChronoSpectra: workspace has no usable bins/processes");

  model.card.SetFlag("workspaces-use-clone", true);
  model.card.ParseDatacard(cfg.datacard, "", "", "", 0, cfg.mass);
  std::set<std::string> card_bins = model.card.bin_set();
  if (card_bins != bins) {
    std::set<std::string> missing = set_difference_names(bins, card_bins);
    std::set<std::string> extra = set_difference_names(card_bins, bins);
    std::ostringstream message;
    message << "ChronoSpectra: workspace/datacard bins are incompatible";
    if (!missing.empty()) message << "; missing in datacard: " << *missing.begin();
    if (!extra.empty()) message << "; extra in datacard: " << *extra.begin();
    throw std::runtime_error(message.str());
  }
  for (std::set<std::string>::const_iterator bin = bins.begin(); bin != bins.end(); ++bin) {
    TH1F reference = model.card.cp().bin(std::vector<std::string>(1, *bin)).GetObservedShape();
    if (reference.GetNbinsX() <= 0) throw std::runtime_error("ChronoSpectra: no usable observed reference shape for bin " + *bin);
    model.references.insert(std::make_pair(*bin, reference));
  }

  model.bin_groups = detail::parse_named_groups(cfg.group_bins, bins, "bin");
  model.process_groups = detail::parse_named_groups(cfg.group_procs, processes, "process");
  validate_name_collisions(cfg, model);
  for (std::map<std::string, detail::NamedGroup>::const_iterator it = model.bin_groups.groups.begin();
       it != model.bin_groups.groups.end(); ++it) {
    model.grouped_bins.insert(it->second.resolved.begin(), it->second.resolved.end());
    TH1F const& reference = model.references.at(it->second.resolved.front());
    for (std::vector<std::string>::const_iterator member = it->second.resolved.begin() + 1;
         member != it->second.resolved.end(); ++member) {
      std::string reason;
      if (!detail::compatible_binning(reference, model.references.at(*member),
                                       it->second.resolved.front(), *member, &reason)) {
        throw std::runtime_error("ChronoSpectra: incompatible bin group " + detail::quote(it->first) +
                                 " between " + detail::quote(it->second.resolved.front()) + " and " +
                                 detail::quote(*member) + ": " + reason);
      }
    }
  }
  for (std::map<std::string, detail::NamedGroup>::const_iterator it = model.process_groups.groups.begin();
       it != model.process_groups.groups.end(); ++it) {
    model.grouped_processes.insert(it->second.resolved.begin(), it->second.resolved.end());
  }

  for (std::set<std::string>::const_iterator bin = bins.begin(); bin != bins.end(); ++bin) {
    if (!model.grouped_bins.count(*bin) || cfg.sep_bin_hists) {
      model.selections.push_back(Selection{*bin, std::vector<std::string>(1, *bin), false, true});
    }
  }
  for (std::map<std::string, detail::NamedGroup>::const_iterator it = model.bin_groups.groups.begin();
       it != model.bin_groups.groups.end(); ++it) {
    model.selections.push_back(Selection{it->first, it->second.resolved, true, false});
  }

  if (cfg.postfit) {
    model.fit.reset(new RooFitResult(ch::OpenFromTFile<RooFitResult>(cfg.fitresult)));
    validate_fit(*model.fit, cfg.fitresult);
  }
  return model;
}

void emit_selection(Config const& cfg, PreparedModel& model, Selection const& selection,
                    detail::OutputManifest& manifest, std::vector<HistogramProduct>& histograms,
                    std::vector<MatrixProduct>& matrices, bool postfit,
                    std::set<std::string>& empty_group_warnings) {
  std::vector<std::string> all_processes = as_vector(select(model.workspace, selection.bins,
                                                            as_vector(model.workspace.process_set())).process_set());
  if (all_processes.empty()) return;
  TH1F const& reference = model.references.at(selection.bins.front());
  std::string prefix = std::string(postfit ? "postfit/" : "prefit/") + selection.name + "/";
  ch::CombineHarvester scope = select(model.workspace, selection.bins, all_processes);

  std::map<std::string, std::vector<std::string> > shape_processes;
  std::vector<std::string> signal_processes = as_vector(scope.cp().signals().process_set());
  std::vector<std::string> background_processes = as_vector(scope.cp().backgrounds().process_set());
  shape_processes["signal"] = signal_processes;
  shape_processes["background"] = background_processes;
  shape_processes["total"] = all_processes;
  for (std::map<std::string, detail::NamedGroup>::const_iterator group = model.process_groups.groups.begin();
       group != model.process_groups.groups.end(); ++group) {
    std::vector<std::string> selected = intersection(all_processes, group->second.resolved);
    if (selected.empty()) {
      std::string warning_key = selection.name + "/" + group->first;
      if (empty_group_warnings.insert(warning_key).second) {
        std::cerr << "ChronoSpectra: warning: process group " << detail::quote(group->first)
                  << " is empty in selection " << detail::quote(selection.name) << "\n";
      }
    } else {
      shape_processes[group->first] = selected;
    }
  }
  for (std::vector<std::string>::const_iterator process = all_processes.begin();
       process != all_processes.end(); ++process) {
    if (!model.grouped_processes.count(*process) || cfg.sep_proc_hists) {
      shape_processes[*process] = std::vector<std::string>(1, *process);
    }
  }

  std::map<std::string, TH1F> restored_shapes;
  for (std::map<std::string, std::vector<std::string> >::const_iterator shape = shape_processes.begin();
       shape != shape_processes.end(); ++shape) {
    if (shape->second.empty()) continue;
    ch::CombineHarvester selected = select(model.workspace, selection.bins, shape->second);
    TH1F evaluated;
    if (postfit) {
      evaluated = cfg.samples == 0 ? selected.GetShapeWithUncertainty()
                                   : selected.GetShapeWithUncertainty(*model.fit, cfg.samples);
    } else {
      evaluated = selected.GetShapeWithUncertainty();
    }
    std::string context = (postfit ? "post-fit " : "pre-fit ") + selection.name + "/" + shape->first;
    restored_shapes[shape->first] = restore_shape(evaluated, reference, context);
    double rate_uncertainty = 0.;
    if (postfit) {
      rate_uncertainty = cfg.samples == 0 ? selected.GetUncertainty()
                                          : selected.GetUncertainty(*model.fit, cfg.samples);
    } else {
      rate_uncertainty = selected.GetUncertainty();
    }
    if (!std::isfinite(rate_uncertainty)) {
      throw std::runtime_error("ChronoSpectra: non-finite total-rate uncertainty in " + context);
    }
    restored_shapes[shape->first].AddBinContent(0,
        rate_uncertainty - restored_shapes[shape->first].GetBinContent(0));
    restored_shapes[shape->first].SetBinError(0, 0.);
    restored_shapes[shape->first].SetTitle(shape->first.c_str());
    add_histogram(manifest, histograms, prefix + shape->first, restored_shapes[shape->first], shape->first);
  }
  ch::CombineHarvester observed_selector = select(model.workspace, selection.bins, all_processes);
  TH1F observed = cfg.skip_obs ? pseudo_observation(restored_shapes.at("total"), cfg.dataset, selection.name)
                               : restore_shape(observed_selector.GetObservedShape(), reference,
                                               (postfit ? "post-fit " : "pre-fit ") + selection.name + "/" + cfg.dataset);
  observed.SetTitle(cfg.dataset.c_str());
  add_histogram(manifest, histograms, prefix + cfg.dataset, observed, cfg.dataset);

  if (postfit && cfg.samples > 0) {
    for (std::map<std::string, std::vector<std::string> >::const_iterator shape = shape_processes.begin();
         shape != shape_processes.end(); ++shape) {
      if (shape->first == cfg.dataset || shape->second.empty()) continue;
      ch::CombineHarvester selected = select(model.workspace, selection.bins, shape->second);
      bool grouped_process = model.grouped_processes.count(shape->first) != 0;
      if (cfg.get_hist_bin_corr && (!grouped_process || cfg.sep_proc_hist_bin_corr) &&
          (!selection.is_actual_bin || !model.grouped_bins.count(selection.name) || cfg.sep_bin_hist_bin_corr)) {
        TH2F correlation = histogram_bin_correlation(selected, *model.fit, cfg.samples,
                                                      selection.name + "/" + shape->first);
        add_matrix(manifest, matrices, prefix + shape->first + "_HistBinCorr", correlation,
                   "Histogram Bin Correlation Matrix");
      }
      if (cfg.get_rate_corr && shape->first != "signal" && shape->first != "background" &&
          shape->first != "total" && model.process_groups.groups.count(shape->first) == 0) continue;
      if (cfg.get_rate_corr && (!selection.is_actual_bin || !model.grouped_bins.count(selection.name) || cfg.sep_bin_rate_corr) &&
          (shape->first == "signal" || shape->first == "background" || shape->first == "total" ||
           model.process_groups.groups.count(shape->first))) {
        TH2F correlation = rate_correlation(selected, *model.fit, cfg.samples,
                                            selection.name + "/" + shape->first);
        add_matrix(manifest, matrices, prefix + shape->first + "_RateCorr", correlation,
                   "Rate Correlation Matrix");
      }
    }
  }
}

void write_products(std::string const& output, detail::OutputManifest const& manifest,
                    std::vector<HistogramProduct>& histograms,
                    std::vector<MatrixProduct>& matrices) {
  (void)manifest;
  fs::path final_path(output);
  fs::path parent = final_path.parent_path();
  if (parent.empty()) parent = fs::path(".");
  if (!fs::exists(parent) || !fs::is_directory(parent)) {
    throw std::runtime_error("ChronoSpectra: output parent directory does not exist: " + parent.string());
  }
  fs::path temporary = parent / fs::unique_path(final_path.filename().string() + ".chrono-%%%%%%.tmp");
  bool committed = false;
  try {
    {
      TFile file(temporary.string().c_str(), "RECREATE");
      if (!file.IsOpen() || file.IsZombie()) throw std::runtime_error("ChronoSpectra: cannot create temporary output " + temporary.string());
      for (std::vector<HistogramProduct>::iterator it = histograms.begin(); it != histograms.end(); ++it) {
        ch::WriteToTFile(&it->histogram, &file, it->path);
      }
      for (std::vector<MatrixProduct>::iterator it = matrices.begin(); it != matrices.end(); ++it) {
        ch::WriteToTFile(&it->matrix, &file, it->path);
      }
      file.Write();
      file.Close();
    }
    TFile validation(temporary.string().c_str(), "READ");
    if (!validation.IsOpen() || validation.IsZombie()) throw std::runtime_error("ChronoSpectra: temporary output failed validation: " + temporary.string());
    validation.Close();
    fs::rename(temporary, final_path);
    committed = true;
  } catch (...) {
    if (!committed && fs::exists(temporary)) fs::remove(temporary);
    throw;
  }
}

int run(Config const& cfg, po::variables_map const& variables) {
  detail::ScopedAddDirectory add_directory(false);
  PreparedModel model = prepare_model(cfg);
  detail::ParameterStateGuard state(model.workspace);
  parse_freeze_tokens(cfg.freeze, model.workspace);

  detail::OutputManifest manifest;
  std::vector<HistogramProduct> histograms;
  std::vector<MatrixProduct> matrices;
  std::vector<SystematicProduct> systematics;
  std::set<std::string> empty_group_warnings;
  if (cfg.store_syst) {
    store_systematics(cfg, model, manifest, histograms, systematics);
  }
  if (!cfg.skipprefit) {
    for (std::vector<Selection>::const_iterator selection = model.selections.begin();
         selection != model.selections.end(); ++selection) {
      emit_selection(cfg, model, *selection, manifest, histograms, matrices, false, empty_group_warnings);
    }
  }
  if (cfg.postfit) {
    model.workspace.UpdateParameters(*model.fit);
    if (variables.count("seed")) {
      unsigned seed = detail::parse_unsigned(variables.at("seed").as<std::string>(), "seed");
      RooRandom::randomGenerator()->SetSeed(seed);
    }
    for (std::vector<Selection>::const_iterator selection = model.selections.begin();
         selection != model.selections.end(); ++selection) {
      emit_selection(cfg, model, *selection, manifest, histograms, matrices, true, empty_group_warnings);
    }
    if (cfg.samples > 0 && cfg.get_rate_corr) {
      TH2F global = rate_correlation(model.workspace, *model.fit, cfg.samples, "global");
      add_matrix(manifest, matrices, "postfit/globalRateCorr", global,
                 "Rate Correlation Matrix");
    }
    add_matrix(manifest, matrices, "postfit/parCorrMat", parameter_correlation(*model.fit),
               "Parameter Correlation Matrix");
  }
  state.restore();
  plot_systematics(cfg, systematics);
  write_products(cfg.output, manifest, histograms, matrices);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    ParsedArguments arguments = parse_arguments(argc, argv);
    if (arguments.help) {
      std::cout << arguments.description << "\n";
      return 0;
    }
    validate_config(arguments.config);
    return run(arguments.config, arguments.variables);
  } catch (po::error const& error) {
    std::cerr << "ChronoSpectra: " << error.what() << "\n";
    return 2;
  } catch (std::exception const& error) {
    std::cerr << "ChronoSpectra: " << error.what() << "\n";
    return 1;
  }
}
