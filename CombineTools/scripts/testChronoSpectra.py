#!/usr/bin/env python3
"""Reproducible end-to-end validation for the SCRAM-built ChronoSpectra tool.

This is intentionally one self-contained Python file. It runs the private
invariant binary, generates a fresh synthetic ROOT shape file and datacard,
converts and fits that fixture, drives the public executable against the
generated workspace, inspects the complete ROOT object tree, exercises
successful and failing CLI paths, checks output transactions and seeded
repeatability, and can write a portable Markdown evidence report. It never
opens analysis-provided datacards, workspaces, or fit results.

The supported invocation is documented in docs/ChronoSpectra.md. The test
expects a CMSSW runtime with PyROOT and uses only the standard Python library
besides ROOT.
"""

import argparse
from array import array
import hashlib
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import time


class TestFailure(AssertionError):
    """An assertion with context suitable for the final test report."""


DEFAULT_FIT = object()
NO_FIT = object()


def check(condition, message):
    if not condition:
        raise TestFailure(message)


def command_text(command):
    return " ".join(str(part) for part in command)


def finite(value):
    return math.isfinite(float(value))


def close(left, right, relative=1e-9, absolute=1e-10):
    return abs(float(left) - float(right)) <= max(
        absolute, relative * max(1.0, abs(float(left)), abs(float(right))))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_root():
    try:
        import ROOT
    except ImportError as error:
        raise TestFailure("PyROOT is required: {}".format(error))
    ROOT.gROOT.SetBatch(True)
    return ROOT


def run_process(command, context, expected=None, cwd=None):
    completed = subprocess.run(
        [str(part) for part in command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=str(cwd) if cwd is not None else None,
    )
    record = {
        "command": command_text(command),
        "returncode": completed.returncode,
        "output": completed.stdout,
        "context": context.current_case or "ChronoSpectra",
    }
    context.commands.append(record)
    if expected is not None and completed.returncode != expected:
        raise TestFailure(
            "{} returned {} (wanted {})\ncommand: {}\n{}".format(
                context.current_case or "ChronoSpectra", completed.returncode,
                expected, command_text(command), completed.stdout))
    return completed


def serialize_root(ROOT, filename):
    """Return a sorted, JSON-safe description of every ROOT key and payload."""
    root_file = ROOT.TFile.Open(str(filename), "READ")
    check(root_file and not root_file.IsZombie(),
          "cannot open ROOT output {}".format(filename))
    entries = []

    def visit(directory, prefix):
        keys = directory.GetListOfKeys()
        for index in range(keys.GetEntries()):
            key = keys.At(index)
            obj = key.ReadObj()
            path = prefix + "/" + key.GetName() if prefix else key.GetName()
            entry = {"path": path, "class": obj.ClassName(), "title": obj.GetTitle()}
            if obj.InheritsFrom("TH1"):
                bins = obj.GetNbinsX()
                entry["bins"] = bins
                entry["edges"] = [obj.GetXaxis().GetBinLowEdge(i)
                                  for i in range(1, bins + 1)]
                if bins:
                    entry["edges"].append(obj.GetXaxis().GetBinUpEdge(bins))
                entry["labels"] = [obj.GetXaxis().GetBinLabel(i)
                                   for i in range(1, bins + 1)]
                entry["content"] = [obj.GetBinContent(i)
                                    for i in range(0, bins + 2)]
                entry["errors"] = [obj.GetBinError(i)
                                   for i in range(0, bins + 2)]
            if obj.InheritsFrom("TH2"):
                xbins = obj.GetNbinsX()
                ybins = obj.GetNbinsY()
                entry["xbins"] = xbins
                entry["ybins"] = ybins
                entry["xlabels"] = [obj.GetXaxis().GetBinLabel(i)
                                    for i in range(1, xbins + 1)]
                entry["ylabels"] = [obj.GetYaxis().GetBinLabel(i)
                                    for i in range(1, ybins + 1)]
                entry["matrix"] = [
                    [obj.GetBinContent(x, y) for y in range(1, ybins + 1)]
                    for x in range(1, xbins + 1)
                ]
            entries.append(entry)
            if obj.InheritsFrom("TDirectory"):
                visit(obj, path)

    try:
        visit(root_file, "")
    finally:
        root_file.Close()
    return sorted(entries, key=lambda item: item["path"])


def entry_map(serialized):
    return {entry["path"]: entry for entry in serialized}


def paths(serialized):
    return set(entry["path"] for entry in serialized)


def require_paths(serialized, required, absent=()):
    available = paths(serialized)
    missing = set(required) - available
    check(not missing, "missing ROOT paths: {}".format(sorted(missing)))
    unexpected = available.intersection(absent)
    check(not unexpected, "unexpected ROOT paths: {}".format(sorted(unexpected)))


def hist(serialized, path):
    entry = entry_map(serialized).get(path)
    check(entry is not None, "missing histogram {}".format(path))
    check(entry.get("class") == "TH1F",
          "{} has class {}, expected TH1F".format(path, entry.get("class")))
    return entry


def matrix(serialized, path):
    entry = entry_map(serialized).get(path)
    check(entry is not None, "missing matrix {}".format(path))
    check(entry.get("class") == "TH2F",
          "{} has class {}, expected TH2F".format(path, entry.get("class")))
    return entry


def compare_json_values(left, right, context):
    if isinstance(left, float) or isinstance(right, float):
        check(close(left, right, relative=1e-11, absolute=1e-12),
              "{} differs: {} != {}".format(context, left, right))
    elif isinstance(left, list):
        check(len(left) == len(right), "{} length differs".format(context))
        for index, (left_value, right_value) in enumerate(zip(left, right)):
            compare_json_values(left_value, right_value,
                                "{}[{}]".format(context, index))
    elif isinstance(left, dict):
        check(left.keys() == right.keys(), "{} keys differ".format(context))
        for key in left:
            compare_json_values(left[key], right[key], "{}/{}".format(context, key))
    else:
        check(left == right, "{} differs: {} != {}".format(context, left, right))


def regular_integral(entry):
    return sum(entry["content"][1:-1])


def meaningfully_different_entries(first, second):
    check(first["bins"] == second["bins"],
          "cannot compare histograms with different bin counts")
    return any(abs(left - right) > 1e-9 * max(1.0, abs(left), abs(right))
               for left, right in zip(first["content"][1:-1],
                                      second["content"][1:-1]))


def check_histogram_contract(serialized, prefix, dataset="data_obs"):
    required = [
        prefix + "/" + dataset,
        prefix + "/signal",
        prefix + "/background",
        prefix + "/total",
    ]
    require_paths(serialized, required)
    objects = entry_map(serialized)
    shape_entries = [objects[path] for path in required]
    for entry in shape_entries:
        check(entry["bins"] > 0, "{} has no regular bins".format(entry["path"]))
        check(len(entry["edges"]) == entry["bins"] + 1,
              "{} has invalid edge count".format(entry["path"]))
        check(len(entry["content"]) == entry["bins"] + 2,
              "{} does not expose both flow cells".format(entry["path"]))
        for value in entry["content"] + entry["errors"]:
            check(finite(value), "{} contains a non-finite value".format(entry["path"]))

    reference = shape_entries[0]
    for entry in shape_entries[1:]:
        check(entry["bins"] == reference["bins"],
              "{} has a different bin count".format(entry["path"]))
        for index, (left, right) in enumerate(zip(entry["edges"], reference["edges"])):
            check(close(left, right, relative=1e-12, absolute=1e-12),
                  "{} edge {} differs".format(entry["path"], index))
        check(entry["labels"] == reference["labels"],
              "{} labels differ".format(entry["path"]))

    signal = objects[prefix + "/signal"]
    background = objects[prefix + "/background"]
    total = objects[prefix + "/total"]
    for index in range(1, total["bins"] + 1):
        check(close(total["content"][index],
                    signal["content"][index] + background["content"][index],
                    relative=1e-6, absolute=1e-5),
              "{} total does not equal signal plus background at bin {}".format(
                  prefix, index))
    for name in ("signal", "background", "total"):
        shape = objects[prefix + "/" + name]
        check(close(shape["errors"][0], 0.0, absolute=1e-12),
              "{} underflow error is not zero".format(shape["path"]))
    check(objects[prefix + "/" + dataset]["class"] == "TH1F",
          "{} observation is not TH1F".format(prefix))
    return reference


def check_matrix_contract(entry, name, require_symmetry=False):
    check(entry["xbins"] == entry["ybins"], "{} is not square".format(name))
    check(len(entry["matrix"]) == entry["xbins"],
          "{} has invalid x dimension".format(name))
    for row in entry["matrix"]:
        check(len(row) == entry["ybins"], "{} has invalid y dimension".format(name))
        for value in row:
            check(finite(value), "{} contains non-finite values".format(name))
            check(-1.000001 <= float(value) <= 1.000001,
                  "{} contains a correlation outside [-1, 1]".format(name))
    if require_symmetry:
        for i in range(entry["xbins"]):
            for j in range(entry["ybins"]):
                check(close(entry["matrix"][i][j], entry["matrix"][j][i],
                            relative=1e-8, absolute=1e-8),
                      "{} is not symmetric at {},{}".format(name, i, j))


def png_dimensions(filename):
    payload = Path(filename).read_bytes()
    check(payload[:8] == b"\x89PNG\r\n\x1a\n" and payload[12:16] == b"IHDR",
          "{} is not a PNG with an IHDR chunk".format(filename))
    return int.from_bytes(payload[16:20], "big"), int.from_bytes(payload[20:24], "big")


def split_fitresult(value):
    check(":" in value, "fit result must have the form file.root:object")
    filename, object_name = value.rsplit(":", 1)
    check(filename and object_name, "fit result has an empty file or object name")
    return filename, object_name


class TestContext:
    def __init__(self, args, root, workspace, datacard, fitresult, work,
                 primary_bin, group_bin_pattern, group_proc_pattern,
                 group_name, group_members, expected_bins, label="synthetic"):
        self.args = args
        self.root = root
        self.workspace = Path(workspace)
        self.datacard = Path(datacard)
        self.fitresult = fitresult
        self.work = Path(work)
        self.cwd = self.datacard.parent.resolve()
        self.primary_bin = primary_bin
        self.group_bin_pattern = group_bin_pattern
        self.group_proc_pattern = group_proc_pattern
        self.group_name = group_name
        self.group_members = list(group_members)
        self.expected_bins = list(expected_bins)
        self.expected_edges = None
        self.expected_labels = None
        self.incompatible_card = None
        self.plot_pattern = getattr(args, "plot_pattern", "catA/sig/shapeN")
        self.label = label
        self.commands = []
        self.cases = []
        self.current_case = ""

    def case(self, name):
        self.current_case = "{}: {}".format(self.label, name)
        self.cases.append(self.current_case)

    def executable_command(self, output, extra=(), workspace=None,
                           datacard=None, fitresult=DEFAULT_FIT):
        command = [
            self.args.executable,
            "--workspace", str(workspace or self.workspace),
            "--datacard", str(datacard or self.datacard),
            "--output", str(output),
        ]
        if fitresult is DEFAULT_FIT:
            selected_fit = self.fitresult if "--postfit" in extra else None
        elif fitresult is NO_FIT:
            selected_fit = None
        else:
            selected_fit = fitresult
        if selected_fit is not None:
            command.extend(["--fitresult", str(selected_fit)])
        command.extend(str(value) for value in extra)
        return command

    def invoke(self, output, extra=(), expected=0, workspace=None,
               datacard=None, fitresult=DEFAULT_FIT):
        return run_process(
            self.executable_command(output, extra, workspace, datacard, fitresult),
            self, expected=expected, cwd=self.cwd)

    def expect_failure(self, name, extra=(), workspace=None,
                       datacard=None, fitresult=DEFAULT_FIT, contains=()):
        self.case(name)
        output = self.work / ("sentinel-{}.root".format(len(self.cases)))
        original = b"ChronoSpectra sentinel: must survive a rejected invocation\n"
        output.write_bytes(original)
        completed = self.invoke(output, extra, expected=None, workspace=workspace,
                                datacard=datacard, fitresult=fitresult)
        check(completed.returncode != 0, "{} unexpectedly succeeded".format(name))
        check("ChronoSpectra:" in completed.stdout,
              "{} emitted no contextual ChronoSpectra diagnostic:\n{}".format(
                  name, completed.stdout))
        check("terminate called" not in completed.stdout and
              "Aborted" not in completed.stdout,
              "{} escaped through an uncaught exception:\n{}".format(
                  name, completed.stdout))
        for text in contains:
            check(text in completed.stdout,
                  "{} diagnostic lacks {!r}:\n{}".format(name, text, completed.stdout))
        check(output.read_bytes() == original,
              "{} replaced the pre-existing output on failure".format(name))
        return completed


def run_internal_invariants(context):
    context.case("private-invariants")
    candidates = []
    if context.args.unit_executable:
        candidates.append(Path(context.args.unit_executable))
    executable = Path(context.args.executable)
    if len(executable.parents) >= 3:
        release = executable.parents[2]
        architecture = executable.parent.name
        candidates.append(release / "test" / architecture / "ChronoSpectraInternalTest")
    candidates.extend([
        executable.with_name("ChronoSpectraInternalTest"),
        Path(os.environ["CHRONOSPECTRA_INTERNAL_TEST"])
        if os.environ.get("CHRONOSPECTRA_INTERNAL_TEST") else Path("/nonexistent"),
    ])
    unit = next((candidate for candidate in candidates if candidate.is_file() and
                 os.access(str(candidate), os.X_OK)), None)
    check(unit is not None,
          "ChronoSpectraInternalTest not found; pass --unit-executable or set "
          "CHRONOSPECTRA_INTERNAL_TEST")
    completed = run_process([str(unit)], context, expected=0)
    check("all invariant tests passed" in completed.stdout,
          "private invariant test did not report completion:\n{}".format(completed.stdout))


def run_help_and_cli_checks(context):
    context.case("help-and-cli-contract")
    completed = run_process([context.args.executable, "--help"], context, expected=0)
    help_text = completed.stdout
    for option in (
            "--workspace", "--datacard", "--output", "--dataset", "--mass",
            "--fitresult", "--postfit", "--skipprefit", "--samples", "--seed",
            "--freeze", "--groupBins", "--groupProcs", "--skipObs",
            "--getRateCorr", "--getHistBinCorr", "--sepProcHists",
            "--sepBinHists", "--sepProcHistBinCorr", "--sepBinHistBinCorr",
            "--sepBinRateCorr", "--storeSyst", "--plotSyst", "--systSaveDir",
            "--logy", "--logLevel", "--log-level"):
        check(option in help_text, "--help does not mention {}".format(option))

    explicit = context.work / "explicit-cli.root"
    context.invoke(explicit, ["--dataset=data_obs", "--mass=125",
                              "--samples=0", "--log-level=warn"], expected=0)
    explicit_tree = serialize_root(context.root, explicit)
    check_histogram_contract(explicit_tree, "prefit/{}".format(context.primary_bin))


def run_parser_edge_checks(context):
    context.case("parser-edge-cases")
    short_help = run_process([context.args.executable, "-h"], context, expected=0)
    check(short_help.stdout == run_process(
        [context.args.executable, "--help"], context, expected=0).stdout,
          "-h and --help produced different help text")

    no_arguments = run_process([context.args.executable], context, expected=None)
    check(no_arguments.returncode != 0, "no arguments unexpectedly succeeded")
    check("ChronoSpectra:" in no_arguments.stdout,
          "no-argument failure lacked a contextual diagnostic")

    missing_output = run_process([
        context.args.executable, "--workspace", str(context.workspace),
        "--datacard", str(context.datacard),
    ], context, expected=None)
    check(missing_output.returncode != 0 and "output" in missing_output.stdout,
          "missing required output option was not rejected")

    short_form = context.work / "short-form.root"
    short_command = [
        context.args.executable, "-w", str(context.workspace),
        "-d", str(context.datacard), "-o", str(short_form),
        "-m", "125", "--samples=0",
    ]
    run_process(short_command, context, expected=0, cwd=context.cwd)
    check_histogram_contract(serialize_root(context.root, short_form),
                             "prefit/{}".format(context.primary_bin))

    for name, extra in (
            ("malformed-boolean", ["--postfit=maybe"]),
            ("negative-samples", ["--samples=-1"]),
            ("nonnumeric-samples", ["--samples=not-a-number"]),
            ("negative-seed", ["--seed=-1"]),
    ):
        context.expect_failure(name, extra, contains=("ChronoSpectra:",))


def run_prefit_checks(context):
    context.case("prefit-shapes-binning-and-observation")
    output = context.work / "prefit.root"
    context.invoke(output, ["--samples=0"], expected=0)
    tree = serialize_root(context.root, output)
    for bin_name in context.expected_bins:
        reference = check_histogram_contract(tree, "prefit/{}".format(bin_name))
        if context.expected_edges is not None:
            check(reference["edges"] == context.expected_edges,
                  "{} does not retain the expected physical edges".format(bin_name))
        if context.expected_labels is not None:
            check(reference["labels"] == context.expected_labels,
                  "{} does not retain the expected bin labels".format(bin_name))
    objects = entry_map(tree)
    for path, entry in objects.items():
        if path.startswith("prefit/"):
            check(entry["class"] in ("TDirectoryFile", "TH1F"),
                  "unexpected pre-fit class at {}: {}".format(path, entry["class"]))

    skip_obs = context.work / "skip-obs.root"
    context.invoke(skip_obs, ["--samples=0", "--skipObs"], expected=0)
    skip_tree = serialize_root(context.root, skip_obs)
    total = hist(skip_tree, "prefit/{}/total".format(context.primary_bin))
    observation = hist(skip_tree, "prefit/{}/data_obs".format(context.primary_bin))
    for index in range(1, total["bins"] + 1):
        check(close(observation["content"][index], total["content"][index]),
              "pseudo-observation differs at regular bin {}".format(index))
        check(finite(observation["errors"][index]) and
              observation["errors"][index] >= 0.,
              "pseudo-observation has invalid Poisson error at bin {}".format(index))
    check(close(observation["content"][0], math.sqrt(regular_integral(total)),
                relative=1e-7, absolute=1e-5),
          "pseudo-observation underflow does not contain sqrt(total yield)")


def run_group_checks(context):
    context.case("groups-and-separation-switches")
    group_options = ["--samples=0",
                     "--groupBins=region:{}".format(context.group_bin_pattern),
                     "--groupProcs={}:{}".format(context.group_name,
                                                  context.group_proc_pattern)]
    grouped = context.work / "grouped.root"
    context.invoke(grouped, group_options, expected=0)
    grouped_tree = serialize_root(context.root, grouped)
    grouped_paths = paths(grouped_tree)
    require_paths(grouped_tree, ["prefit/region/total",
                                 "prefit/region/{}".format(context.group_name)])
    check("prefit/{}/total".format(context.primary_bin) not in grouped_paths,
          "grouped actual-bin output was not suppressed")
    for member in context.group_members:
        check("prefit/region/{}".format(member) not in grouped_paths,
              "grouped member process {} was not suppressed".format(member))

    proc_separated = context.work / "separate-processes.root"
    context.invoke(proc_separated, group_options + ["--sepProcHists"], expected=0)
    proc_tree = serialize_root(context.root, proc_separated)
    require_paths(proc_tree, ["prefit/region/{}".format(member)
                              for member in context.group_members])
    check("prefit/{}/total".format(context.primary_bin) not in paths(proc_tree),
          "--sepProcHists unexpectedly restored grouped actual-bin output")

    bin_separated = context.work / "separate-bins.root"
    context.invoke(bin_separated, group_options + ["--sepBinHists"], expected=0)
    bin_tree = serialize_root(context.root, bin_separated)
    require_paths(bin_tree, ["prefit/{}/total".format(context.primary_bin),
                             "prefit/region/{}".format(context.group_name)])
    check("prefit/region/{}".format(context.group_members[0]) not in paths(bin_tree),
          "--sepBinHists unexpectedly restored grouped process members")

    all_separated = context.work / "all-separated.root"
    context.invoke(all_separated, group_options + ["--sepProcHists", "--sepBinHists"],
                   expected=0)
    all_tree = serialize_root(context.root, all_separated)
    require_paths(all_tree, ["prefit/{}/total".format(context.primary_bin),
                             "prefit/{}/{}".format(context.primary_bin,
                                                    context.group_members[0]),
                             "prefit/region/total",
                             "prefit/region/{}".format(context.group_name)])

    malformed = [
        ("skipprefit-without-postfit", ["--skipprefit"], DEFAULT_FIT),
        ("postfit-without-fitresult", ["--postfit"], NO_FIT),
        ("store-syst-without-prefit", ["--postfit", "--skipprefit", "--storeSyst"], DEFAULT_FIT),
        ("plot-without-store-syst", ["--plotSyst=all"], DEFAULT_FIT),
        ("invalid-log-level", ["--logLevel=verbose"], DEFAULT_FIT),
        ("unknown-option", ["--does-not-exist"], DEFAULT_FIT),
        ("invalid-group-regex", ["--groupBins=bad:["], DEFAULT_FIT),
        ("unmatched-group-regex", ["--groupBins=bad:not_a_real_bin"], DEFAULT_FIT),
        ("unsafe-group-name", ["--groupBins=unsafe/name:{}".format(context.primary_bin)], DEFAULT_FIT),
        ("duplicate-group-name", ["--groupBins=x:{};x:{}".format(
            context.primary_bin, context.primary_bin)], DEFAULT_FIT),
        ("missing-dataset", ["--dataset=not_in_workspace"], DEFAULT_FIT),
    ]
    for name, extra, fitresult in malformed:
        context.expect_failure(name, extra, fitresult=fitresult,
                               contains=("ChronoSpectra:",))
    if context.incompatible_card is not None:
        context.expect_failure(
            "incompatible-group-binning",
            ["--samples=0", "--groupBins=region:cat.*"],
            datacard=context.incompatible_card,
            contains=("catA", "catB"))


def fit_parameter_names(ROOT, fitresult):
    filename, object_name = split_fitresult(fitresult)
    root_file = ROOT.TFile.Open(filename, "READ")
    check(root_file and not root_file.IsZombie(),
          "cannot open fit result {}".format(filename))
    fit = root_file.Get(object_name)
    check(fit, "fit result object {} is missing".format(fitresult))
    parameters = fit.floatParsFinal()
    names = [parameters.at(index).GetName() for index in range(parameters.getSize())]
    root_file.Close()
    check(names, "fit result has no floating parameters")
    return names


def run_postfit_checks(context):
    context.case("postfit-unsampled-and-state-options")
    postfit_zero = context.work / "postfit-zero.root"
    context.invoke(postfit_zero, ["--postfit", "--samples=0"],
                   fitresult=context.fitresult, expected=0)
    zero_tree = serialize_root(context.root, postfit_zero)
    for bin_name in context.expected_bins:
        check_histogram_contract(zero_tree, "prefit/{}".format(bin_name))
        check_histogram_contract(zero_tree, "postfit/{}".format(bin_name))
    require_paths(zero_tree, ["postfit/parCorrMat"])
    zero_paths = paths(zero_tree)
    check(not any(path.endswith("RateCorr") or path.endswith("HistBinCorr")
                  for path in zero_paths),
          "samples=0 emitted a sampling-derived matrix")
    check(not close(regular_integral(hist(zero_tree,
                                         "prefit/{}/total".format(context.primary_bin))),
                    regular_integral(hist(zero_tree,
                                          "postfit/{}/total".format(context.primary_bin))),
                    absolute=1e-6),
          "combined post-fit output did not retain distinct pre-fit and post-fit values")

    par = matrix(zero_tree, "postfit/parCorrMat")
    check(par["xlabels"] and par["xlabels"] == par["ylabels"],
          "parameter-correlation axes are not labelled consistently")
    check_matrix_contract(par, "postfit/parCorrMat", require_symmetry=True)

    skipped = context.work / "postfit-skipped-prefit.root"
    context.invoke(skipped, ["--postfit", "--skipprefit", "--samples=0"],
                   fitresult=context.fitresult, expected=0)
    skipped_paths = paths(serialize_root(context.root, skipped))
    check(not any(path.startswith("prefit/") for path in skipped_paths),
          "--skipprefit emitted pre-fit products")
    check("postfit/{}/total".format(context.primary_bin) in skipped_paths,
          "--skipprefit omitted post-fit products")

    parameter = fit_parameter_names(context.root, context.fitresult)[0]
    frozen = context.work / "freeze.root"
    context.invoke(frozen, ["--samples=0", "--mass=125",
                            "--freeze={}".format(parameter)], expected=0)
    assigned = context.work / "freeze-assigned.root"
    context.invoke(assigned, ["--samples=0", "--freeze={} = 0".format(parameter).replace(" ", "")],
                   expected=0)

    missing_fit = "{}:not_a_fit_result".format(split_fitresult(context.fitresult)[0])
    context.expect_failure("missing-fit-object", ["--postfit"], fitresult=missing_fit,
                           contains=("not_a_fit_result",))
    context.expect_failure("missing-workspace-file", [],
                           workspace=context.work / "does-not-exist.root",
                           contains=("does-not-exist.root",))
    context.expect_failure("missing-datacard-file", [],
                           datacard=context.work / "does-not-exist.txt",
                           contains=("does-not-exist.txt",))
    empty_workspace = context.work / "empty-workspace.root"
    empty_file = context.root.TFile.Open(str(empty_workspace), "RECREATE")
    empty_file.Close()
    context.expect_failure("missing-workspace-object", [], workspace=empty_workspace,
                           contains=("object 'w'",))


def run_sampling_checks(context):
    context.case("sampling-correlation-switches-and-seed")
    sampled = context.work / "sampled.root"
    sample_options = ["--postfit", "--samples=1", "--seed=17"]
    context.invoke(sampled, sample_options, fitresult=context.fitresult, expected=0)
    sampled_tree = serialize_root(context.root, sampled)
    sampled_paths = paths(sampled_tree)
    require_paths(sampled_tree, ["postfit/globalRateCorr",
                                 "postfit/{}/total_RateCorr".format(context.primary_bin),
                                 "postfit/{}/total_HistBinCorr".format(context.primary_bin)])
    for path in sorted(sampled_paths):
        if path.endswith("RateCorr"):
            check_matrix_contract(matrix(sampled_tree, path), path)
        if path.endswith("HistBinCorr"):
            check_matrix_contract(matrix(sampled_tree, path), path, require_symmetry=True)

    repeated = context.work / "sampled-repeat.root"
    context.invoke(repeated, sample_options, fitresult=context.fitresult, expected=0)
    compare_json_values(sampled_tree, serialize_root(context.root, repeated),
                        "fixed-seed output")

    disabled = context.work / "correlations-disabled.root"
    context.invoke(disabled, ["--postfit", "--samples=1", "--getRateCorr=false",
                              "--getHistBinCorr=false"], fitresult=context.fitresult,
                   expected=0)
    disabled_paths = paths(serialize_root(context.root, disabled))
    check(not any(path.endswith("RateCorr") or path.endswith("HistBinCorr")
                  for path in disabled_paths),
          "correlation switches did not suppress all matrices")

    rate_disabled = context.work / "rate-correlation-disabled.root"
    context.invoke(rate_disabled, ["--postfit", "--samples=1",
                                   "--getRateCorr=false"],
                   fitresult=context.fitresult, expected=0)
    rate_disabled_tree = serialize_root(context.root, rate_disabled)
    rate_disabled_paths = paths(rate_disabled_tree)
    check("postfit/globalRateCorr" not in rate_disabled_paths and
          not any(path.endswith("RateCorr") for path in rate_disabled_paths),
          "getRateCorr=false did not suppress rate matrices")
    check(any(path.endswith("HistBinCorr") for path in rate_disabled_paths),
          "getRateCorr=false incorrectly suppressed histogram-bin matrices")

    hist_disabled = context.work / "hist-correlation-disabled.root"
    context.invoke(hist_disabled, ["--postfit", "--samples=1",
                                   "--getHistBinCorr=false"],
                   fitresult=context.fitresult, expected=0)
    hist_disabled_tree = serialize_root(context.root, hist_disabled)
    hist_disabled_paths = paths(hist_disabled_tree)
    check("postfit/globalRateCorr" in hist_disabled_paths and
          any(path.endswith("RateCorr") for path in hist_disabled_paths),
          "getHistBinCorr=false incorrectly suppressed rate matrices")
    check(not any(path.endswith("HistBinCorr") for path in hist_disabled_paths),
          "getHistBinCorr=false did not suppress histogram-bin matrices")

    group_options = [
        "--postfit", "--samples=1", "--seed=17",
        "--groupBins=region:{}".format(context.group_bin_pattern),
        "--groupProcs={}:{}".format(context.group_name, context.group_proc_pattern),
        "--sepProcHists", "--sepBinHists", "--sepProcHistBinCorr",
        "--sepBinHistBinCorr", "--sepBinRateCorr",
    ]
    separated = context.work / "grouped-correlations.root"
    context.invoke(separated, group_options, fitresult=context.fitresult, expected=0)
    separated_tree = serialize_root(context.root, separated)
    require_paths(separated_tree, [
        "postfit/region/{}_HistBinCorr".format(context.group_name),
        "postfit/region/{}_RateCorr".format(context.group_name),
        "postfit/{}/{}_HistBinCorr".format(context.primary_bin, context.group_name),
        "postfit/{}/{}_RateCorr".format(context.primary_bin, context.group_name),
    ])


def run_systematic_checks(context):
    context.case("systematics-and-plot-transaction")
    output = context.work / "systematics.root"
    plot_directory = context.work / "plots"
    context.invoke(output, ["--samples=0", "--storeSyst",
                            "--plotSyst={}".format(context.plot_pattern),
                            "--systSaveDir", str(plot_directory), "--logy"], expected=0)
    tree = serialize_root(context.root, output)
    syst_paths = [path for path in paths(tree) if path.startswith("systematics/")]
    check(any(path.endswith("_Up") for path in syst_paths),
          "--storeSyst emitted no Up variation")
    check(any(path.endswith("_Down") for path in syst_paths),
          "--storeSyst emitted no Down variation")
    if context.label == "synthetic":
        nominal = hist(tree, "systematics/catA/sig")
        pure_up = hist(tree, "systematics/catA/sig_syst/shapeN_Up")
        pure_down = hist(tree, "systematics/catA/sig_syst/shapeN_Down")
        check(close(regular_integral(nominal), regular_integral(pure_up),
                    relative=1e-6, absolute=1e-6),
              "synthetic pure-shape Up variation changed its regular integral")
        check(close(regular_integral(nominal), regular_integral(pure_down),
                    relative=1e-6, absolute=1e-6),
              "synthetic pure-shape Down variation changed its regular integral")
        check(meaningfully_different_entries(nominal, pure_up) and
              meaningfully_different_entries(nominal, pure_down),
              "synthetic pure-shape variations were not retained")
    plots = sorted(plot_directory.glob("*.png"))
    check(plots, "--plotSyst emitted no PNG files")
    for plot in plots:
        width, height = png_dimensions(plot)
        check(width >= 2500 and height >= 2200,
              "systematic plot is not high resolution: {}x{}".format(width, height))

    unmatched = context.work / "unmatched-pattern.root"
    unmatched_result = context.invoke(
        unmatched, ["--samples=0", "--storeSyst", "--plotSyst=no/no/no",
                     "--systSaveDir", str(context.work / "no-plots")], expected=0)
    check("matched no stored variation" in unmatched_result.stdout,
          "unmatched systematic pattern did not produce a warning")

    existing_plot = plot_directory / plots[0].name
    original_plot = existing_plot.read_bytes()
    context.expect_failure(
        "plot-collision", ["--samples=0", "--storeSyst",
                            "--plotSyst={}".format(context.plot_pattern),
                            "--systSaveDir", str(plot_directory)])
    check(existing_plot.read_bytes() == original_plot,
          "plot collision changed an existing plot")


def run_io_failure_checks(context):
    context.case("output-transaction")
    output_parent = context.work / "missing-parent"
    check(not output_parent.exists(), "test setup unexpectedly contains output parent")
    output = output_parent / "output.root"
    completed = run_process(
        context.executable_command(output, ["--samples=0"]), context, expected=None,
        cwd=context.cwd)
    check(completed.returncode != 0, "missing output parent unexpectedly succeeded")
    check("output parent directory" in completed.stdout,
          "missing output parent diagnostic is not contextual")
    check(not output.exists(), "ChronoSpectra created output under missing parent")


def _write_shape(ROOT, root_file, name, values, labels,
                 bin_edges=(0.0, 0.5, 2.0, 5.0)):
    edges = array("d", bin_edges)
    shape = ROOT.TH1F(name, name, len(bin_edges) - 1, edges)
    shape.SetDirectory(root_file)
    for index, value in enumerate(values, 1):
        shape.SetBinContent(index, float(value))
        shape.SetBinError(index, math.sqrt(abs(float(value))))
        shape.GetXaxis().SetBinLabel(index, labels[index - 1])
    shape.Write()


def generate_synthetic_fixture(context):
    """Create a compact two-channel shape card and fit result in a private dir."""
    context.case("synthetic-fixture-generation")
    fixture = context.work / "synthetic-fixture"
    fixture.mkdir(parents=True, exist_ok=True)
    labels = ["low", "middle", "high"]
    shapes_path = fixture / "synthetic_shapes.root"
    shape_values = {
        "catA": {
            "sig": [4.0, 3.0, 2.0],
            "bkgA": [8.0, 7.0, 6.0],
            "bkgB": [3.0, 4.0, 5.0],
            "data_obs": [15.0, 14.0, 13.0],
        },
        "catB": {
            "sig": [3.0, 4.0, 5.0],
            "bkgA": [9.0, 6.0, 5.0],
            "bkgB": [4.0, 5.0, 6.0],
            "data_obs": [16.0, 15.0, 16.0],
        },
    }
    root_file = context.root.TFile.Open(str(shapes_path), "RECREATE")
    check(root_file and not root_file.IsZombie(), "could not create synthetic shape file")
    try:
        for category, processes in shape_values.items():
            for process, values in processes.items():
                _write_shape(context.root, root_file,
                             "{}_{}".format(category, process), values, labels)
            for process, values in processes.items():
                if process == "data_obs":
                    continue
                total = sum(values)
                if process == "sig":
                    up = [4.4, 3.0, 1.6] if category == "catA" else [3.3, 4.0, 4.7]
                    down = [3.6, 3.0, 2.4] if category == "catA" else [2.7, 4.0, 5.3]
                    check(close(sum(up), total) and close(sum(down), total),
                          "synthetic pure-shape variation changed its integral")
                elif process == "bkgA":
                    up = [1.10 * value for value in values]
                    down = [0.90 * value for value in values]
                else:
                    up = [value + 0.25 for value in values]
                    down = [max(0.05, value - 0.25) for value in values]
                _write_shape(context.root, root_file,
                             "{}_{}_shapeNUp".format(category, process), up, labels)
                _write_shape(context.root, root_file,
                             "{}_{}_shapeNDown".format(category, process), down, labels)
        root_file.Write()
    finally:
        root_file.Close()

    card = fixture / "synthetic_card.txt"
    card.write_text("""imax 2 number of channels
jmax 2 number of backgrounds
kmax * number of nuisance parameters
------------------------------------------------------------
shapes * catA synthetic_shapes.root catA_$PROCESS catA_$PROCESS_$SYSTEMATIC
shapes * catB synthetic_shapes.root catB_$PROCESS catB_$PROCESS_$SYSTEMATIC
------------------------------------------------------------
bin catA catB
observation -1 -1
------------------------------------------------------------
bin catA catA catA catB catB catB
process sig bkgA bkgB sig bkgA bkgB
process 0 1 2 0 1 2
rate 9 21 12 12 20 15
------------------------------------------------------------
lumi lnN 1.05 1.05 1.05 1.05 1.05 1.05
bkgNorm lnN - 1.10 1.10 - 1.10 1.10
shapeN shape 1 1 1 1 1 1
""", encoding="utf-8")
    workspace = fixture / "synthetic_workspace.root"
    check(context.args.text2workspace,
          "text2workspace.py is required for synthetic tests")
    run_process([context.args.text2workspace, str(card), "-m", "125", "-o", str(workspace)],
                context, expected=0, cwd=fixture)
    check(workspace.is_file(), "text2workspace did not create synthetic workspace")

    check(context.args.combine, "combine is required to create synthetic fit results")
    run_process([context.args.combine, "-M", "FitDiagnostics", str(workspace), "-m", "125",
                 "-n", ".chronospectra_synthetic", "--saveWorkspace",
                 "--cminDefaultMinimizerStrategy=0"], context, expected=0, cwd=fixture)
    fit_path = fixture / "fitDiagnostics.chronospectra_synthetic.root"
    check(fit_path.is_file(), "combine did not create synthetic fit result")
    fit_file = context.root.TFile.Open(str(fit_path), "READ")
    check(fit_file and not fit_file.IsZombie(), "cannot open synthetic fit result")
    try:
        fit = fit_file.Get("fit_s")
        check(fit and fit.InheritsFrom("RooFitResult"),
              "synthetic fit result does not contain fit_s")
    finally:
        fit_file.Close()
    incompatible_shapes = fixture / "incompatible_shapes.root"
    incompatible_file = context.root.TFile.Open(str(incompatible_shapes), "RECREATE")
    check(incompatible_file and not incompatible_file.IsZombie(),
          "could not create incompatible synthetic shape file")
    try:
        for category, processes in shape_values.items():
            edges = (0.0, 1.0, 2.0, 5.0) if category == "catB" else (0.0, 0.5, 2.0, 5.0)
            for process, values in processes.items():
                _write_shape(context.root, incompatible_file,
                             "{}_{}".format(category, process), values, labels,
                             bin_edges=edges)
            for process, values in processes.items():
                if process == "data_obs":
                    continue
                if process == "sig":
                    up = [4.4, 3.0, 1.6] if category == "catA" else [3.3, 4.0, 4.7]
                    down = [3.6, 3.0, 2.4] if category == "catA" else [2.7, 4.0, 5.3]
                elif process == "bkgA":
                    up = [1.10 * value for value in values]
                    down = [0.90 * value for value in values]
                else:
                    up = [value + 0.25 for value in values]
                    down = [max(0.05, value - 0.25) for value in values]
                _write_shape(context.root, incompatible_file,
                             "{}_{}_shapeNUp".format(category, process), up, labels,
                             bin_edges=edges)
                _write_shape(context.root, incompatible_file,
                             "{}_{}_shapeNDown".format(category, process), down, labels,
                             bin_edges=edges)
        incompatible_file.Write()
    finally:
        incompatible_file.Close()
    incompatible_card = fixture / "incompatible_card.txt"
    incompatible_card.write_text(card.read_text(encoding="utf-8").replace(
        "synthetic_shapes.root", "incompatible_shapes.root"), encoding="utf-8")
    return fixture, workspace, card, "{}:fit_s".format(fit_path), incompatible_card


def run_public_suite(context):
    run_help_and_cli_checks(context)
    run_parser_edge_checks(context)
    run_prefit_checks(context)
    run_group_checks(context)
    run_postfit_checks(context)
    run_sampling_checks(context)
    run_systematic_checks(context)
    run_io_failure_checks(context)


def write_report(context, report_path, started, completed):
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    source_root = Path(__file__).resolve()
    target_executable = Path(context.args.executable).resolve()
    repository = source_root.parents[2]
    try:
        commit = subprocess.run(
            ["git", "-C", str(repository), "rev-parse", "HEAD"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.strip()
    except OSError:
        commit = "unavailable"
    command_parts = [
        "python3", str(source_root),
        "--executable", str(context.args.executable),
    ]
    if context.args.unit_executable:
        command_parts.extend(["--unit-executable", str(context.args.unit_executable)])
    if context.args.skip_internal:
        command_parts.append("--skip-internal")
    if context.args.text2workspace:
        command_parts.extend(["--text2workspace", str(context.args.text2workspace)])
    if context.args.combine:
        command_parts.extend(["--combine", str(context.args.combine)])
    command_parts.extend(["--report", str(report_path)])
    reproduction = " ".join(shlex.quote(str(part)) for part in command_parts)
    lines = [
        "# ChronoSpectra test execution report",
        "",
        "Status: **{}**".format("PASS" if completed else "FAIL"),
        "",
        "The single harness CombineTools/scripts/testChronoSpectra.py generated "
        "all fixture inputs in a temporary directory and ran the public "
        "integration checks against those generated inputs{}. No analysis "
        "datacard, workspace, or fit result was opened, and the source checkout "
        "was not modified.".format(
            " and the private invariant binary" if not context.args.skip_internal else ""),
        "",
        "## Reproduction",
        "",
        "    {}".format(reproduction),
        "",
        "## Environment",
        "",
        "| Field | Value |",
        "| --- | --- |",
        "| executable | {} |".format(target_executable),
        "| ChronoSpectra SHA-256 | {} |".format(sha256(target_executable)),
        "| harness SHA-256 | {} |".format(sha256(source_root)),
        "| repository commit | {} |".format(commit),
        "| Python | {} |".format(sys.version.split()[0]),
        "| platform | {} |".format(platform.platform()),
        "| generated inputs | temporary fixture (removed after the run) |",
        "| started | {} |".format(time.ctime(started)),
        "| duration | {:.2f}s |".format(time.time() - started),
        "",
        "## Coverage executed",
        "",
        "- private regex/path/glob, online-correlation, binning/flow, global ROOT-state, and parameter-restoration invariants;" if not context.args.skip_internal else
        "- private invariant binary not run in this fallback invocation (use the SCRAM test gate);",
        "- help, aliases, required options, explicit mass/dataset/log-level, and all public option declarations;",
        "- pre-fit/post-fit and post-fit-only products, physical binning, labels, flow cells, aggregate arithmetic, pseudo-data, and fit-parameter correlations;",
        "- regex bin/process groups, duplicate/unmatched/unsafe/colliding definitions, overlap suppression, and all five separation switches;",
        "- sampled rate/histogram-bin correlations, matrix shape/labels/finite-range/symmetry contracts, fixed-seed structural repeatability, and disabled-correlation switches;",
        "- freeze forms, missing input/object/fit cases, invalid combinations, contextual diagnostics, uncaught-exception protection, and sentinel preservation;",
        "- systematic storage, variation discovery, wildcard no-match warnings, high-resolution/log-scale PNGs, plot collision safety, and output-parent failures.",
        "",
        "## Command results",
        "",
        "| Case | Exit |",
        "| --- | ---: |",
    ]
    for command in context.commands:
        lines.append("| {} | {} |".format(command["context"], command["returncode"]))
    lines.extend([
        "",
        "Cases executed: **{}**. Commands executed: **{}**.".format(
            len(context.cases), len(context.commands)),
        "",
        "Command output is retained in memory only to avoid placing generated logs "
        "or ROOT products in the source repository. The full command line, input "
        "paths, executable hash, source hash, and environment are recorded here.",
        "",
    ])
    lines.extend([
        "The fixture-generation commands include PyROOT shape creation, "
        "text2workspace.py conversion, combine FitDiagnostics, and construction "
        "of an incompatible-binning negative-test card. All generated files were "
        "temporary.",
        "",
    ])
    report_path.write_text("\n".join(lines), encoding="utf-8")


def default_report_path():
    requested = os.environ.get("CHRONOSPECTRA_REPORT")
    if requested:
        return requested
    candidate = (Path(__file__).resolve().parents[3] /
                 "chronospectra_assimilation_planning" /
                 "CHRONOSPECTRA_GENERATED_FIXTURE_REPORT.md")
    return str(candidate) if candidate.parent.is_dir() else None


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", default=shutil.which("ChronoSpectra"),
                        help="ChronoSpectra executable; defaults to PATH")
    parser.add_argument("--unit-executable",
                        help="ChronoSpectraInternalTest path; defaults beside --executable")
    parser.add_argument("--skip-internal", action="store_true",
                        help="skip the private invariant binary (diagnostic fallback only)")
    parser.add_argument("--text2workspace",
                        default=shutil.which("text2workspace.py"),
                        help="text2workspace.py used to build the generated workspace")
    parser.add_argument("--combine", default=shutil.which("combine"),
                        help="combine executable used to build the generated fit result")
    parser.add_argument("--report", default=None,
                        help="write a Markdown execution report; defaults to planning output when available")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    started = time.time()
    root = load_root()
    if args.executable is None:
        raise TestFailure("ChronoSpectra executable not found on PATH; pass --executable")
    if args.report is None:
        args.report = default_report_path()
    work = Path(tempfile.mkdtemp(prefix="chronospectra-validation-"))
    generation = TestContext(
        args, root, work / "placeholder.root", work / "placeholder.txt",
        "placeholder.root:fit_s", work, "catA", "cat.*", "bkg.*", "bkg",
        ["bkgA", "bkgB"], ["catA", "catB"], label="fixture-generation")
    context = generation
    completed = False
    try:
        fixture, workspace, card, fitresult, incompatible_card = generate_synthetic_fixture(
            generation)
        context = TestContext(
            args, root, workspace, card, fitresult, work / "synthetic-suite",
            "catA", "cat.*", "bkg.*", "bkg", ["bkgA", "bkgB"],
            ["catA", "catB"], label="synthetic")
        context.expected_edges = [0.0, 0.5, 2.0, 5.0]
        context.expected_labels = ["low", "middle", "high"]
        context.incompatible_card = incompatible_card
        context.work.mkdir(parents=True, exist_ok=True)
        context.commands = generation.commands
        context.cases = generation.cases
        if not args.skip_internal:
            run_internal_invariants(context)
        run_public_suite(context)
        completed = True
        if args.report:
            write_report(context, args.report, started, completed)
        print("testChronoSpectra.py: all {} cases and {} commands passed".format(
            len(context.cases), len(context.commands)))
        return 0
    except (TestFailure, AssertionError) as error:
        if args.report:
            try:
                write_report(context, args.report, started, completed)
            except Exception as report_error:
                print("testChronoSpectra.py: could not write failure report: {}".format(
                    report_error), file=sys.stderr)
        print("testChronoSpectra.py: FAIL: {}".format(error), file=sys.stderr)
        return 1
    finally:
        if os.environ.get("CHRONOSPECTRA_KEEP_TEST_WORK"):
            print("testChronoSpectra.py: retained temporary work at {}".format(work))
        else:
            shutil.rmtree(str(work), ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
