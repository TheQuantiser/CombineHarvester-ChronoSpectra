ChronoSpectra {#chronospectra}
==============================

ChronoSpectra is a SCRAM-built CombineHarvester executable that evaluates a
Combine workspace and its original datacard and writes structured one-
dimensional pre-fit and post-fit shape products. It composes the existing
CombineHarvester workspace, datacard, uncertainty, correlation, and binning
APIs. The existing PostFitShapesFromWorkspace and
PostFit2DShapesFromWorkspace tools are unchanged.

The executable adds:

* named bin and process groups with full-match regular expressions;
* restored physical bin edges and labels, including compatible grouped bins;
* pre-fit and post-fit histograms with target-defined uncertainties;
* sampled rate and local histogram-bin correlations;
* fit-parameter correlation output;
* systematic variation storage and PNG plots; and
* contextual diagnostics with transactional ROOT and plot output.

Prerequisites
-------------

Use the CMSSW and CombinedLimit versions supported by the current
CombineHarvester release. The validated baseline is:

| Component | Version |
| --- | --- |
| CMSSW | CMSSW_16_0_0 |
| SCRAM architecture | el9_amd64_gcc13 |
| CombinedLimit | v11.0.0 |
| ROOT | 6.36.07 |
| GCC | 13.4.0 |

The workspace must contain:

* a RooWorkspace named w;
* a RooStats ModelConfig whose PDF is a RooSimultaneous;
* the selected dataset, data_obs by default; and
* at least one usable bin and process.

The datacard is required because it supplies the observed reference shapes
used to restore regular-bin edges and labels.

Usage
-----

~~~text
ChronoSpectra --workspace workspace.root --datacard card.txt \
  --output chronospectra.root
~~~

The complete option list is emitted by the executable itself:

~~~bash
ChronoSpectra --help
~~~

| Option | Default | Meaning |
| --- | --- | --- |
| --workspace, -w | required | ROOT file containing w |
| --datacard, -d | required | Original datacard for mass and binning |
| --output, -o | required | ROOT output, committed only after success |
| --dataset | data_obs | Workspace dataset and observation key |
| --mass, -m | 125 | Mass passed to ParseDatacard |
| --fitresult, -f | empty | ROOT file and object: file.root:object |
| --postfit | false | Add post-fit products |
| --skipprefit | false | Omit pre-fit products; requires --postfit |
| --samples | 2000 | Post-fit samples; zero uses target unsampled errors |
| --seed | unset | RooFit seed set before the first sampling operation |
| --freeze | empty | Comma-separated PARAM or PARAM=VALUE entries |
| --groupBins | empty | Named full-match regular-expression bin groups |
| --groupProcs | empty | Named full-match regular-expression process groups |
| --skipObs | false | Emit pseudo-data from the total prediction |
| --getRateCorr | true | Emit sampled rate matrices |
| --getHistBinCorr | true | Emit sampled local histogram-bin matrices |
| --sepProcHists | false | Restore grouped individual-process histograms |
| --sepBinHists | false | Restore grouped actual-bin histograms |
| --sepProcHistBinCorr | false | Restore grouped process bin matrices |
| --sepBinHistBinCorr | false | Restore grouped actual-bin bin matrices |
| --sepBinRateCorr | false | Restore grouped actual-bin rate matrices |
| --storeSyst | false | Store pre-fit individual-process variations |
| --plotSyst | empty | all or comma-separated bin/process/parameter patterns |
| --systSaveDir | shapeSystPlots | PNG destination |
| --logy | false | Use a logarithmic upper plot panel when valid |
| --logLevel, --log-level | info | info, warn, or error |

Boolean options accept a bare spelling and an explicit =true or =false
spelling. Unknown options, malformed values, unsafe names, invalid
combinations, missing ROOT objects, invalid fit results, and output
preflight/I/O errors return nonzero with a ChronoSpectra: diagnostic.
--help returns zero before required file checks.

Operating modes
---------------

| skipprefit | postfit | samples | Result |
| ---: | ---: | ---: | --- |
| false | false | any | Pre-fit products |
| false | true | 0 | Pre-fit, unsampled post-fit, and parCorrMat |
| false | true | >0 | Pre-fit and sampled post-fit products |
| true | true | 0 | Unsampled post-fit products and parCorrMat |
| true | true | >0 | Sampled post-fit products |
| true | false | any | Invalid |

Post-fit always requires a valid RooFitResult, including samples=0. Sampling-
derived rate and histogram-bin matrices are omitted at zero samples.
parCorrMat is emitted for every successful post-fit run.

Groups
------

Groups use the following grammar:

~~~text
name:member,member;other_name:member
~~~

Names and members are trimmed. Each member is a Boost regular expression
matched against the complete model name. Use .*fragment.* for substring
matching. Duplicate members are normalized. Invalid expressions, empty
tokens, unmatched expressions, unsafe names, duplicate names, and collisions
with model names, aggregate keys, the dataset, or another group are rejected
before output creation.

Overlapping groups are independent sums. Individual outputs covered by a
group are suppressed by default and can be restored with the corresponding
sep* option. Actual names and group names are traversed lexicographically.
A grouped bin is valid only when all member reference shapes have the same
regular-bin count, labels, and edges within the documented floating-point
tolerance.

Output products
---------------

For each selected actual bin or compatible bin group, the standard keys are:

~~~text
prefit/<selection>/data_obs
prefit/<selection>/signal
prefit/<selection>/background
prefit/<selection>/total
prefit/<selection>/<process-group>
prefit/<selection>/<process>

postfit/<selection>/data_obs
postfit/<selection>/signal
postfit/<selection>/background
postfit/<selection>/total
postfit/<selection>/<process-group>
postfit/<selection>/<process>
~~~

All shape products are TH1F. signal and background are omitted when empty;
total is required for a nonempty selection. The observation key is always
present. With skipObs=false it is the observed shape. With skipObs=true it is
pseudo-data copied from total, with ROOT Poisson errors and a historical
sqrt(total regular yield) summary in its underflow cell. Negative or
non-finite pseudo-data predictions are rejected.

Regular-bin edges and labels are restored from the datacard reference shape.
Evaluated underflow and overflow cells are preserved. For model shapes, the
underflow content is the target CombineHarvester total-rate uncertainty and
its error is zero; ROOT Integral() continues to exclude flow cells. Negative
regular model bins are not globally clipped or replaced by absolute values.

Statistics and matrices
-----------------------

Unsampled uncertainties use the current target definitions. Sampled shape-bin
errors use population variance around the sample mean while retaining the
nominal central content. Sampled total-rate uncertainty remains the target
nominal-centered quantity.

When postfit and samples>0 are selected, eligible products include:

~~~text
postfit/<selection>/<key>_RateCorr
postfit/<selection>/<key>_HistBinCorr
postfit/globalRateCorr
~~~

Rate matrices come from the unchanged target rate-correlation API. Only
non-finite cells at the output boundary are replaced by zero, with one
contextual warning per affected matrix. Finite target values are not
rewritten.

Local histogram-bin correlations use the sample-mean-centered population
covariance:

~~~text
Cij = (1/N) sum_k (hki - mean_i) (hkj - mean_j)
rhoij = Cij / sqrt(Cii Cjj)
~~~

The implementation uses an online Welford/Chan accumulator and retains
O(B^2) state rather than all sampled histograms. Zero-variance cells are zero,
positive-variance diagonal cells are one, and the matrix is symmetric.
Disable this calculation with --getHistBinCorr=false when the quadratic work
is not needed.

parCorrMat is copied from RooFitResult::correlationMatrix() in
floatParsFinal() order and is present for every post-fit run.

Systematics and plots
---------------------

storeSyst requires pre-fit output. Stored products are:

~~~text
systematics/<actual-bin>/<process>
systematics/<actual-bin>/<process>_syst/<parameter>_Up
systematics/<actual-bin>/<process>_syst/<parameter>_Down
~~~

Frozen and bin-wise no-effect parameters are skipped. Pure-shape changes are
retained even when their regular-bin integral is unchanged. plotSyst selects
variations with exactly three glob components: bin/process/parameter. all
and */*/* select everything; .* is accepted as a compatibility spelling for
*. Valid unmatched patterns warn and produce no plot. Plots are PNG only,
rendered in a staging directory, checked for collisions, and moved to the
requested systSaveDir after successful rendering.

With logy enabled, the upper panel is logarithmic only when all plotted
regular-bin values are positive. Otherwise the upper panel remains linear
with a warning.

Output and state safety
-----------------------

ChronoSpectra takes one outer snapshot of parameter values, asymmetric
errors, ranges, frozen flags, groups, and linked RooFit variables. Requested
freezes, fit updates, systematic shifts, and sampling mutations are restored
on success and during exception unwinding. TH1::AddDirectory is also restored.

The complete output manifest is built before the requested output is touched.
Products are written to a unique temporary ROOT file beside the destination,
reopened for validation, closed, and then renamed into place. A rejected or
failed invocation leaves an existing output untouched and removes only its
own temporary file. Plot staging follows the same owned-file principle.

Testing
-------

The reproducible test driver is the single Python file

~~~text
CombineTools/scripts/testChronoSpectra.py
~~~

It uses PyROOT and the standard Python library. The normal invocation creates
every test input itself in a temporary directory: a ROOT shape file, a text
datacard, a `text2workspace.py` workspace, a `combine` FitDiagnostics result,
and an incompatible-binning card for a negative test. It never discovers or
opens an analysis-provided datacard, workspace, fit result, or input ROOT file.
The generated fixture has nonuniform edges, bin labels, compatible grouped
bins, two backgrounds, a rate nuisance, and both rate-changing and pure-shape
variations. No generated ROOT file or datacard is checked into the repository.
The driver automatically finds the SCRAM-installed ChronoSpectraInternalTest
beside the CMSSW release test directory when the production executable was
built in the standard layout. Pass `--unit-executable` explicitly when using
another layout.
`--skip-internal` exists only for diagnosing a runtime that cannot launch the
SCRAM test binary; it is not a replacement for the private invariant gate.

The test has two layers:

1. ChronoSpectraInternalTest directly checks the private regex, path, glob,
   online-correlation, binning/flow, TH1 global-state, and parameter-
   restoration invariants, including exception unwinding.
2. testChronoSpectra.py drives fresh subprocesses and inspects ROOT products
   and PNGs. It checks help and aliases; pre-fit, post-fit, post-fit-only,
   samples=0, and seeded sampled modes; physical binning and flow cells;
   aggregate shapes and pseudo-data; freeze forms; bin/process groups;
   all five sep* switches; rate and histogram-bin matrix contracts; matrix
   suppression switches; systematics and plotting; invalid input/object/fit
   paths; contextual diagnostics; sentinel preservation; plot collisions;
   missing output parents; and fixed-seed structural repeatability.

Run inside the supported CMSSW runtime from `$CMSSW_BASE/src`:

~~~bash
export SCRAM_ARCH=el9_amd64_gcc13
eval "$(scram runtime -sh)"
scram b -j 8
scram b runtests

python3 CombineHarvester/CombineTools/scripts/testChronoSpectra.py \
  --executable "$(command -v ChronoSpectra)" \
  --unit-executable "$CMSSW_BASE/test/$SCRAM_ARCH/ChronoSpectraInternalTest"
~~~

That command generates its own inputs, converts and fits them, runs the private
invariant binary, and runs the full public suite against the generated
workspace and fit result. The generated fixture is created, converted, fitted,
tested, and removed under the harness temporary directory. Use
`--text2workspace` or `--combine` only when those tools are not found on
`PATH`; use `--report` or `CHRONOSPECTRA_REPORT` to choose the report location.

The official CMS EL9 route for a host that cannot run the release directly is:

~~~bash
/cvmfs/cms.cern.ch/common/cmssw-el9 --command-to-run '
  set -e
  export SCRAM_ARCH=el9_amd64_gcc13
  cd "$CMSSW_BASE/src"
  eval "$(scram runtime -sh)"
  scram b -j 8
  scram b runtests
  python3 CombineHarvester/CombineTools/scripts/testChronoSpectra.py \
    --executable "$(command -v ChronoSpectra)" \
    --unit-executable "$CMSSW_BASE/test/$SCRAM_ARCH/ChronoSpectraInternalTest"
'
~~~

The test creates a temporary directory for products, plots, and generated
fixture inputs. Set
CHRONOSPECTRA_KEEP_TEST_WORK=1 to retain that directory for debugging. The
report records the executable and harness hashes, repository commit, runtime,
fixture-generation commands, command count, case count, and pass/fail result.
In the checked-out workspace it is written automatically below the planning
directory; set `CHRONOSPECTRA_REPORT` or pass `--report` to choose another
location. Generated ROOT products should remain outside the Git worktree.

The report records the commands used to construct the generated workspace and
fit result, and the executable and harness hashes. These results are evidence
for the exact runtime used; they are not a claim of cross-ROOT bitwise
reproducibility.

Coverage boundary
-----------------

The single-file runner provides black-box coverage for the implemented
executable against a reproducible generated multi-channel fixture. The
generated campaign covers multiple compatible channels, physical binning and
labels, pure-shape and rate-changing variations, workspace conversion, a
fitted post-fit result, grouping, sampling, incompatible grouped-bin
rejection, and output transactions. Negative regular-bin inputs, zero-variance
and mixed-variance models, autoMCStats, wrong-class ModelConfig/PDF objects,
and full large-B/N performance sweeps remain separate targeted campaigns.
The official 2D extractor must be regression-tested separately because
ChronoSpectra itself produces one-dimensional outputs. These boundaries are
reported explicitly rather than inferred as passing from a nearby smoke test.

Examples
--------

The following are direct ChronoSpectra CLI examples. They are not inputs to
the synthetic test driver; that driver generates its own workspace, datacard,
and fit result as described above.

~~~bash
# Pre-fit products with datacard binning
ChronoSpectra -w workspace.root -d card.txt -o shapes.root

# Sampled post-fit products with a reproducible stream
ChronoSpectra -w workspace.root -d card.txt -o shapes.root \
  -f fitDiagnostics.root:fit_s --postfit --samples=2000 --seed=17

# Named groups and pseudo-data
ChronoSpectra -w workspace.root -d card.txt -o grouped.root \
  --groupBins='control:.*control.*;signal:signal_region' \
  --groupProcs='diboson:.*ZZ.*' --skipObs --samples=0
~~~

Known limitations and attribution
----------------------------------

RooFit sampling uses a process-global random generator. A seed provides
repeatability only for the same executable, ROOT/CMSSW environment, input,
seed, sample count, and complete option set. It is not a promise of bitwise
stability across ROOT releases or different traversal options. Dense local
bin matrices are quadratic in the number of regular bins.

ChronoSpectra concept and original implementation: Mohammad Abrar Wadud
(2024). This CombineHarvester integration is a clean-room composition of the
target public APIs. Exact attribution and license wording for upstream
submission remains subject to maintainer and original-author review.
