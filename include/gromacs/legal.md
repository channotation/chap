## Vendored GROMACS Headers ##

CHAP currently vendors selected GROMACS-derived headers to allow its compilation with newer GROMACS versions (>=2018) that have changed their API. This includes are mostly reimplementations of random number generators of other libraries with additional upstream notices:

1. `random/threefry.h` (includes notice for the Boost Software License)
2. `random/uniformrealdistribution.h` (includes notice for The LLVM compiler infrastructure MIT-licensed implementation)

To stay compliant when vendoring these files, use the following checklist:

1. Place vendored files in a dedicated folder, for example: `/third_party/gromacs/`.
2. Do not modify these files unless you are prepared to distribute those modifications under the GNU Lesser General Public License (LGPL) v2.1 or later, consistent with GROMACS licensing terms.
3. In your About/Legal material, include this exact statement: "This software uses the ThreeFry implementation from GROMACS, licensed under LGPL v2.1."
4. Provide a link to the GROMACS source code: https://github.com/gromacs/gromacs

When redistributing, preserve all original license and attribution notices present in the vendored files, including the additional Boost and LLVM-related notices where applicable.