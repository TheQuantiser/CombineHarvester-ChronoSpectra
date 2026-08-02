# Standalone CombineHarvester build

The standalone build is additive. The existing CMSSW/SCRAM files and targets
remain authoritative and are not modified by CMake. The supported standalone
dependency environment is the preinstalled `combine-11` conda environment.

The build requires an explicit staging prefix:

```bash
conda run -n combine-11 env -u LD_LIBRARY_PATH cmake \
  -S /home/mohammad/CH/CombineHarvester \
  -B /tmp/chronospectra-dual-build \
  -G Ninja \
  -DCH_CONDA_PREFIX=/home/mohammad/.conda/envs/combine-11 \
  -DHiggsAnalysisCombinedLimit_ROOT=/home/mohammad/.conda/envs/combine-11 \
  -DCMAKE_INSTALL_PREFIX=/tmp/chronospectra-dual-install
conda run -n combine-11 env -u LD_LIBRARY_PATH cmake \
  --build /tmp/chronospectra-dual-build --parallel 8
conda run -n combine-11 env -u LD_LIBRARY_PATH ctest \
  --test-dir /tmp/chronospectra-dual-build --output-on-failure
conda run -n combine-11 env -u LD_LIBRARY_PATH cmake \
  --install /tmp/chronospectra-dual-build
```

CombinedLimit is never fetched, rebuilt, copied, or added as a subdirectory.
`FindHiggsAnalysisCombinedLimit.cmake` validates its header, shared library,
PCM, and rootmap under the supplied prefix. ROOT, Boost, Eigen3, LibXml2, and
VDT are also re-found from the approved environment.

The build also requires `patchelf` on `PATH`. The conda-forge GCC driver used
by `combine-11` embeds that environment's absolute library path in every link;
the CMake post-link step replaces it on the three CH artifacts with the
relative paths shown below, so the installed package can be used from a
CMSSW-selected runtime.

The standalone targets are:

* `CombineHarvester::CombineTools`, the complete current `CombineTools/src`
  library;
* `ChronoSpectra`, the existing executable compiled without source changes;
* `ChronoSpectraInternalTest`, the existing focused invariant test.

The ROOT dictionary is built from the unchanged
`CombineTools/src/classes.h` and `classes_def.xml` inputs. ROOT 6.36
recursively emits a collection proxy for `ch::CombineHarvester`'s private
`AutoMCStatsSettings` map even though the class is transient, which makes the
generated C++ fail access checking. CMake therefore writes a build-tree-only
ROOT metadata header using `ROOT::Meta::Selection::kNoAutoSelected` for
`auto_stats_settings_`. The outer `ch::CombineHarvester` dictionary remains
selected, while ROOT does not auto-select that private member type. This is
dictionary metadata, not a replacement C++ declaration: the library sources,
public/private API, CMSSW selection file, and numerical implementation are
unchanged.

Because `rootcling --noIncludePaths` intentionally omits build/source paths
from the dictionary metadata, ROOT runtime checks should expose the installed
include root explicitly, for example:

```bash
ROOT_INCLUDE_PATH=/tmp/chronospectra-dual-install/include
LD_LIBRARY_PATH=/tmp/chronospectra-dual-install/lib:$CONDA_PREFIX/lib
ROOT_LIBRARY_PATH=/tmp/chronospectra-dual-install/lib:$CONDA_PREFIX/lib
```

The installed executable has an `$ORIGIN/../lib` RUNPATH for the CH library;
the explicit `LD_LIBRARY_PATH` above supplies the external ROOT, Boost, VDT,
and CombinedLimit runtime libraries without depending on an inherited
historical loader path.

The staged install contains the library, executable, public headers, CMake
package metadata, and uniquely named ROOT dictionary artifacts:

```text
lib/libCombineHarvesterCombineTools.so
lib/libCombineHarvesterCombineTools_rdict.pcm
lib/libCombineHarvesterCombineTools.rootmap
bin/ChronoSpectra
include/CombineHarvester/CombineTools/interface/
lib/cmake/CombineHarvester/
```

To consume an installed tree, keep the same dependency environment active and
point `CMAKE_PREFIX_PATH` at the staging prefix. The package config re-finds
the external dependencies and does not embed the source or build directory.
