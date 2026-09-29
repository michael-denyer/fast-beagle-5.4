# Notice

fast-beagle is licensed under GPL-3.0-or-later ([LICENSE](LICENSE)). It is a derivative work of Beagle 5.5 and 5.4 and contains code from the projects below. Each entry names the files, the upstream copyright and the upstream license.

## Beagle 5.4 and 5.5

- Files: the Java source in `java/src/`, Beagle 5.4 (29Oct24), and every file under `src/` ported from Beagle 5.5 (27Feb25) or 5.4.
- Copyright (C) 2014-2024 Brian L. Browning.
- License: GPL-3.0-or-later.

`java/src/` holds the Beagle 5.4 Java source unmodified. Each file under `src/` that is ported from a Beagle source file carries Browning's copyright line and a notice that names the Java file and the year of modification.

## PLINK 2.0

- Files: `src/bgen/bgen_writer.c`, `src/bgen/bgen_writer.h`, `src/bgen/plink2_num.c` and `src/bgen/plink2_num.h`.
- Version: PLINK 2.0 v2.0.0-a.7.8.
- Copyright (C) 2005-2026 Shaun Purcell, Christopher Chang.
- License: GPL-3.0-or-later for `plink2_import.cc`, `plink2_data.cc`, `plink2_filter.cc` and `plink2_export.cc`. LGPL-3.0-or-later for `plink2_pvar.cc` and `include/plink2_string.cc`.

`src/bgen/bgen_writer.*` ports code from `plink2_import.cc`, `plink2_data.cc`, `plink2_pvar.cc`, `plink2_filter.cc` and `plink2_export.cc`. `src/bgen/plink2_num.*` ports code from `include/plink2_string.cc`. LGPL-3.0-or-later permits conveying this code as part of a GPL-3.0-or-later work. Each of the four files carries the PLINK copyright line, names the upstream files, and states that it was modified in 2026.

## fdlibm 5.3

- Files: `src/jcompat/fdlibm/e_log.c`, `e_log10.c`, `e_pow.c` and `s_expm1.c`.
- Source: [netlib fdlibm 5.3](https://www.netlib.org/fdlibm/).
- Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
- License: a permissive notice. Each file keeps it.

The notice reads: "Permission to use, copy, modify, and distribute this software is freely granted, provided that this notice is preserved." Each file also states its 2026 modifications. `src/jcompat/fdlibm/words.h` is new code, not part of fdlibm.

## libdeflate 1.25

- Files: `third_party/libdeflate/`.
- Source: plink2's source tree at tag v2.0.0-a.7.8.
- Copyright 2016 Eric Biggers. Copyright 2024 Google LLC.
- License: MIT. See [third_party/libdeflate/COPYING](third_party/libdeflate/COPYING).

## htslib

fast-beagle links htslib at build time. htslib distributes its files outside `cram/` under the MIT/Expat license, and its `cram/` files under a modified 3-clause BSD license.

`src/main/vcf_index.c` ports the VCF branch of `tbx_parse1` from htslib's `tbx.c`, which htslib does not export.

- Copyright (C) 2009, 2010, 2012-2015, 2017-2020, 2022-2023, 2025-2026 Genome Research Ltd.
- Copyright (C) 2010-2012 Broad Institute.
- License: MIT. The license text in `tbx.c` follows.

```text
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
```
