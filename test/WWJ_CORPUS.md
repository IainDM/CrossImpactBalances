# The Weimer-Jehle corpus

38 ScenarioWizard models and 21 of ScenarioWizard's own solution sets, supplied by
**Wolfgang Weimer-Jehle** — the author of the Cross-Impact Balance method and of ScenarioWizard —
from his own working files, in August 2026.

They are the strongest validation this package has. 19 of the models ship with the consistent
scenarios ScenarioWizard itself found; on every one of those, JuCIB's independent exhaustive search
returns the **identical set** — not a matching count, the same scenarios — across spaces from 10¹⁰ to
10²⁴. The remaining 19 models have no published solution, and JuCIB's kernels for them are recorded
here as its own result.

## Redistribution

**These files are not in this repository.** They are Weimer-Jehle's, and no permission to publish
them has been given. `test/wwj_corpus/` is listed in `.gitignore`.

> **Licence / redistribution terms:** _to be filled in once agreed with Wolfgang Weimer-Jehle._

Everything *derived* from the corpus is committed: the test code (`test/wwj_tests.jl`), the verifier
(`test/verify_wwj.jl`), and the machine-readable results (`test/bench_results_wwj.json`). The tests
detect that the data is absent and skip, so a clone without it is green.

## How it is exercised

| | What runs | Cost |
|---|---|---|
| `Pkg.test()` | All 38 parsed and shape-checked; all 19 ScenarioWizard solution sets confirmed to be genuine fixed points; **17 of the 19 pairs re-searched and compared set-for-set**; the stale exports and blank-variant behaviour pinned | ~1m55s, 1 thread |
| `JUCIB_WWJ_FULL=1 Pkg.test()` | Adds N40a and D55b (the two slowest pairs), the 19 models with no reference solution, and the composed 10²⁷/10³⁰ kernels | ~7 min, 1 thread |
| `julia -t auto --project=. test/verify_wwj.jl` | Everything, both disagreement directions, non-zero exit on any discrepancy. **The acceptance test.** | ~25 s, 8 threads |

The gap between the first two rows is `--check-bounds=yes`, which `Pkg.test()` sets and which disables
every `@inbounds` in the search kernels — about 7× on these searches. N40 takes 4.8 s in an ordinary
session and 34.8 s under the test runner; N40a 10.7 s becomes 72.8 s and D55b 14.3 s becomes 99.2 s.
Those last two alone are 59% of the cost of re-searching all 19, which is why they are the two held
back. ScenarioWizard's answers *for* them are still checked by default — only JuCIB's independent
re-search of them is deferred.

### Staging a copy

Put the 38 `.scw` and 21 `.sl` files directly in `test/wwj_corpus/` — flat, no subdirectories, keeping
their original names and bytes. Then:

```bash
julia -t auto --project=. test/verify_wwj.jl
```

Check the files against the manifest at the end of this document first; the results below were
measured against exactly those bytes.

### When redistribution is agreed

1. Add `test/wwj_corpus/** -text` to `.gitattributes` **before** staging the files for commit. Without
   it, `core.autocrlf` normalises 8 MB of CRLF in the index, the committed bytes stop being the bytes
   Weimer-Jehle sent, and the manifest below passes on Windows and fails on Linux. Do not add a rule
   for `test/sample_files` at the same time — that would rewrite existing fixture bytes.
2. Delete the `test/wwj_corpus/` line from `.gitignore` and `git add` the corpus (8.17 MB raw,
   about +0.36 MB packed).
3. Re-run `Pkg.test()`. `test/wwj_tests.jl` stops skipping with no code change.
4. Land the `paper/paper.md` wording, which is deliberately held back while a reviewer cannot obtain
   the inputs, and record the terms above.

## What is in the files

`.scw` files carry `$ ScenarioWizard 4.0`; the `.sl` files carry `$ ScenarioWizard 3.0`. That version
mismatch is in the bytes as supplied and is not a problem — the solution-list format did not change.

The models are **anonymised**. Descriptors are named `A.`, `B.` … `Z.`, then `a` … `x`, then `51` …
`100`; variants are `A1`, `A2`, `A3` and so on. Past the printable-ASCII run the name generator walks
off the end of the alphabet into `[1`, `\1`, `]1`, `^1` and then produces **blank names** — the bare
` -` lines in the variant blocks are that artefact, not missing data. Every file has some; N20 has 4,
B100 has 140. Nothing of the underlying studies is present, and no file has a duplicate descriptor
name.

Blank variant names are worth knowing about because they make *name-based* addressing ambiguous:
`get_impact(cib, "e", "", …)` resolves to the first match rather than raising. Index-based addressing
is unaffected, and that is what everything here uses. `test/wwj_tests.jl` pins the behaviour so it
cannot change silently.

## The corpus

Kernel sizes and search times are JuCIB's, single machine, `julia -t auto` (i7-9700, 8 threads),
2026-08-10, via `test/verify_wwj.jl`. Whole corpus: **24.6 s**.

| Model | Desc | Var | Scenarios | Kernel | ScenarioWizard `.sl` | Search |
|---|---:|---:|---:|---:|---|---|
| N20 | 20 | 54 | 10^8.4 | 33 | none | exhaustive |
| N25 | 25 | 67 | 10^10.4 | 3 | **3** — identical set | exhaustive |
| N25a | 25 | 67 | 10^10.4 | 3 | none | exhaustive |
| N25b | 25 | 67 | 10^10.4 | 2 | none | exhaustive |
| N25c | 25 | 67 | 10^10.4 | 10 | none | exhaustive |
| D30a | 30 | 75 | 10^11.7 | 13 | **13** — identical set | exhaustive |
| N30 | 30 | 81 | 10^12.6 | 15 | **15** — identical set | exhaustive |
| D35 | 35 | 88 | 10^13.7 | 51 | **51** — identical set | exhaustive |
| D35a | 35 | 88 | 10^13.7 | 15 | none | exhaustive |
| N35 | 35 | 94 | 10^14.6 | 3 | **3** — identical set | exhaustive |
| D40 | 40 | 100 | 10^15.6 | 110 | **110** — identical set | exhaustive |
| D40a | 40 | 100 | 10^15.6 | 3,240 | none | split_cib |
| D40_a | 40 | 120 | 10^18.6 | 379 | none | exhaustive |
| N40 | 40 | 108 | 10^16.8 | 15 | **15** — identical set | exhaustive |
| N40a | 40 | 108 | 10^16.8 | 18 | **18** — identical set | exhaustive |
| D45 | 45 | 113 | 10^17.6 | 86 | **86** — identical set | exhaustive |
| D45a | 45 | 113 | 10^17.6 | 12 | none | exhaustive |
| N45 | 45 | 121 | 10^18.8 | 62 | **stale** (29) — see below | exhaustive |
| N45a | 45 | 121 | 10^18.8 | 2 | none | exhaustive |
| N45b | 45 | 121 | 10^18.8 | 18 | none | exhaustive |
| N45c | 45 | 121 | 10^18.8 | 42 | none | exhaustive |
| B50 | 50 | 100 | 10^15.1 | 172 | **172** — identical set | exhaustive |
| D50 | 50 | 125 | 10^19.5 | 470 | **470** — identical set | exhaustive |
| D50a | 50 | 125 | 10^19.5 | 35 | none | exhaustive |
| N50 | 50 | 135 | 10^21.0 | 11 | **stale** (93) — see below | exhaustive |
| D55 | 55 | 138 | 10^21.5 | 6,340 | none | exhaustive |
| D55a | 55 | 138 | 10^21.5 | 140 | **140** — identical set | exhaustive |
| D55b | 55 | 138 | 10^21.5 | 328 | **328** — identical set | exhaustive |
| D55c | 55 | 138 | 10^21.5 | 620 | **620** — identical set | exhaustive |
| B60 | 60 | 120 | 10^18.1 | 928 | **928** — identical set | exhaustive |
| D60 | 60 | 150 | 10^23.3 | 194 | **194** — identical set | exhaustive |
| D60b | 60 | 120 | 10^18.1 | 382 | **382** — identical set | exhaustive |
| B70 | 70 | 140 | 10^21.1 | 3,156 | **3,156** — identical set | exhaustive |
| B80 | 80 | 160 | 10^24.1 | 18,432 | **18,432** — identical set | exhaustive |
| B90 | 90 | 180 | 10^27.1 | 18,874,368 | none | split_cib |
| B100 | 100 | 200 | 10^30.1 | 19,327,352,832 | none | split_cib |
| B100b | 100 | 200 | 10^30.1 | 1,195,648 | none | split_cib |
| B100c | 100 | 200 | 10^30.1 | 26,992 | none | exhaustive |

**19 paired models, 25,136 ScenarioWizard scenarios, zero disagreements** — nothing ScenarioWizard
found that JuCIB missed, and nothing JuCIB found that ScenarioWizard missed.

**14 models exceed `typemax(Int64)` scenarios** and so have no signatures at all: B70, B80, B90, B100,
B100b, B100c, D50, D50a, D55, D55a, D55b, D55c, D60, N50. `find_consistent` searches them anyway —
branch-and-bound never numbers a scenario — while `max_signature` and `find_basins` refuse, as they
must. This corpus is why that distinction exists.

**Four models have kernels too large to hold as a list.** B90's is 1.9×10⁷ scenarios and B100's is
1.9×10¹⁰; a direct search for B90's passed 15 GB of resident memory before being killed. Their
influence maps split into independent islands, so `split_cib` computes the size exactly as a product
and never builds it — B100 in 1.9 s. `test/verify_wwj.jl` takes that route whenever a model has more
than one island and nothing needs the scenarios themselves.

## Known-bad files

Three of the 22 `.sl` files Weimer-Jehle supplied match no model in the corpus. This is recorded so
nobody re-derives it.

**`B90.sl` — not staged, not committed.** 22 descriptors wide with variant indices up to 4, so it is
not a binary model and cannot be `B90.scw` (90 binary descriptors); no 22-descriptor model exists in
the corpus at all. Its header declares 1,241,136 solutions but it holds 100,000 — a truncated export.
At 5.90 MB it is 42% of the raw corpus and 66% of its compressed cost, for no test value, so it is
excluded from `test/wwj_corpus/` entirely. `test/wwj_tests.jl` asserts it is absent.

**`N45.sl` and `N50.sl` — staged, and pinned as failing.** Both have the right width and index range
for `N45.scw` / `N50.scw`, so they load cleanly, but **none** of their scenarios is a fixed point of
that matrix: 0 of 29 for N45, 0 of 93 for N50. JuCIB's kernels for those models are 62 and 11. The
straightforward reading is that they were exported from an earlier revision of the matrices, which
were then edited. They are kept because they are useful in both directions — a change that made them
start passing would be as suspicious as one that made a good pair start failing — and both counts are
pinned in `test/wwj_tests.jl`.

## Manifest

SHA-256 of every staged file, against which the results above were measured.

```
5132ebb9f3039f0d32c1398d23ed5566dad67f92f122d637659ce40d4b13c09f   312626  B100.scw
0b9edf65a5954147f012d28c6cf8a0ed930109977ed36ab775fdde0aa5ebb7ca   308306  B100b.scw
63b730950c4ca7882c6d159a25be46a93929e1122d03855603522699be58743e   308306  B100c.scw
9f694bd86f95b24a4250c49c7f46bf446e3628f64c2dd73989cfb74d7067261e    77778  B50.scw
8357e0d22fef5bb1e9b2f5a50d9394baf0a2dc916d8229195a6855e8c1f1c85f    19822  B50.sl
37d2f47c1b8385dd1e37df34d2200ce69e716322cacfea9de4ce9ae23c734145   112498  B60.scw
48640b1498b8e6c86ee331e259c2bb74f667be6458a632640d99d3fbd2bb464b   125322  B60.sl
6fc60c82883eb75059986dd253bf47f4ebc3a91267481d3d38cd38204a44ecdb   153618  B70.scw
73e9bbe7087705bdbb81782c9adaac2e548344fae7eceeba32fe107cfefdad67   489222  B70.sl
051bba090165f5348c4def3e2d9e6a95c036c6b7934fcf0a6569d5dece97075b   201138  B80.scw
8c8959202b5f5b4f52922e73adc87e1499e14f6b52d921be7571c3ecd8682dfc  3225643  B80.sl
7c8d8d50521279577b1eb757ad92e7008ea686c2e7283fed9b50020571dfbdfa   254384  B90.scw
2e1030f6b20150159f15f29bb4b9fba88bed500554e2bbf2956035ac352e4f64    37202  D30a.scw
38a3e7644cd24015232d116964d07ad9d0a874983111f136542316d9f98a18dd     1017  D30a.sl
13a4c8d2052e74db6b5cc24473283d1733aee1c56b1d4964865a305297d45664    49758  D35.scw
385aad551179b7ec1e5b37c15018a587d6c652e0a817dbe6225f73107ad951dc     4377  D35.sl
502390534e25d0cb43ba4899ddc4c5e4c8659c8f3fbfaa528ea3bd74cbdc175c    50157  D35a.scw
1d497ebd4f760d01adb389a80d1d1381da5d68b91b035d10c8a0c007813feb90    64248  D40.scw
80fb61070fe06a7b580d47cd8a3efb0600b58264f8eca89fa83429205643c571    10492  D40.sl
02a00fa5386e2f9ae4f4da23210b680d41cec050695929ae493a6e8fb5762be2    82168  D40_a.scw
5bd2431bb454d77599ecae5350bdb9434f1c3876323a05e48d6c3d3935983cd8    65118  D40a.scw
da59e8adbe14f905029bccf02698d5de0ca108c764955ac9f29c587d6ca7a326    81007  D45.scw
5171ef2212a33816821b7daec37fc890f2c0f23bf1f515b6b1e16d0a75136283     9072  D45.sl
ea532284fe48bcd8f3634a31495c448d5b17b5310a14002b34803ca3a3a87247    81007  D45a.scw
0b2136c1909a35063d7fd664feb6e81242ee78ee1dcd137a94d3df77cb059a40    98698  D50.scw
5e41805d6f84cae66c68cf4069d457b89f344304a74368ac0d1d55dff5d4aa2c    54092  D50.sl
588dd434318f4d3f2a8acf211de7a13a2c9a67d98583a7f4809e97e0bbbd066d    98698  D50a.scw
c5a7b45ed5dc91df17b20574baf8e6237233a78e745bfb5fea80fd6fb7c3ea9e   122188  D55.scw
704a8e41786eb1677db3fa2708587f28f4872c22e70b364c9016e83490bdfe49   119262  D55a.scw
f4b14ccb32a8145cd67cb69ae0e187aa8f4613c473e59d0a3da6de3a6eeb35a2    17542  D55a.sl
f5a66943e0aa45e3fa059d985e41ba5250563fbc5690e8f0d7036aad76014e18   119262  D55b.scw
9055f489b170e5244a8a9ca17b6e825d20c43c86e42638e206550daf551a26d1    41042  D55b.sl
b148b294c87f83e93ac01381fed0a415e9596edba2062607234241905ae05f0b   119262  D55c.scw
b39cf3a8597b25b60d58c2e6d33cf30895f5372748bc47c77efa37a19db332ca    77542  D55c.sl
02d67c98933ea1f08897341f63cefb9baf6c2a11bedaff30eaf6d41ade655599   140558  D60.scw
a074eafee05c5d92fcf602b1df44a104fc679ab240d20947edddfc9a01ddad6f    26232  D60.sl
9309ea2a316df9baa8985279e66849559dbf1810a3d9dcab60b251af7d004d45   107618  D60b.scw
2995d10bb569f9821864004b6211448738690ab149d4efd2c36f8217a22fd156    51612  D60b.sl
0f5039059c1166ce15f52d0da07004995df4e4b3f056cc198f018dbb21b73590    18758  N20.scw
63dd07cbe4dfcd0b2f448bcb9450ea8f3179af7d6ab68006d6cc49b8aa23d984    28542  N25.scw
6c5db89ca1f732ea3d755be0c01e8d2043b0f7e2793e2e35107ad5179186ce87      237  N25.sl
0b7f094eb7ba0cf542db05dbc33c190637c5ceb04746a1d54e2c3e8387b84170    28542  N25a.scw
479a05561579deb4bcf5e3b940048ffa27f4218672d2fe9cbf266f4415a2fecd    28542  N25b.scw
01775eb073b437139c0cf4c418f62a1d254b4bcbb20fc395558dc5dfc1c64f02    28542  N25c.scw
555d8ea330d9fa4a8cf5a218f1ddb3668d60fd20dfa546ef89c881bc79a396b6    40684  N30.scw
922ce605c8e6f9043494a428b23d9cddc4a4bec83f20fd817cdb9f7e81664856     1167  N30.sl
5a6187e7fa6f1ed697c3851d98c119381312c25d970ee2eddf6163b267d4746a    54222  N35.scw
6076f5470ecded7ee1075a770d8bf4d823cea08a2294510ee90c0e2ac4b141a6      297  N35.sl
027302b75487ab57b1e7374e0256a71de9f72a36898edbf7ac22d1814e133b96    70588  N40.scw
af4d75de53ea35829bf2dd2537d15a73c3d29e2e47ffa4f4275e7c10740a45fc     1467  N40.sl
301e2d61761cf33a678ad6fcf38bc079189dd658df47918db01aba0f55ff72be    70588  N40a.scw
cbac734c017ce72bbae0316427c18a566a4c44ded9bf206a3993fa1e7e0878db     1752  N40a.sl
a25a2b913bae68fee1767abef9976aa5fc95cd35ade655c63c4c2bf6a13b45cc    88134  N45.scw
3174e30f54d028cdcefa6f32005c3bccd30bedae32079662cb3d0ae9cadb3097     3087  N45.sl
e448506526389181e5d968defdb87ee36218ba138a2bdde16f8e539e6eb46b87    88134  N45a.scw
6111903fa3260ae74b169b18b519a337f534b2779771c54f6fa2ac71e17cab93    88134  N45b.scw
4d028f95e0675faf7041e8f76247979e756f6a70953e1325cb3b6e64e3a7470f    88134  N45c.scw
30ca64eb3aa4f41cbae405a095518d0b0ee8111f22070020b70faedeaa140e3f   109258  N50.scw
58ba3415aab9a9c8884a957cd783f64bdad9c8a6184a5cf6609db4ea55da6eda    10731  N50.sl
```

`B90.sl` is deliberately not listed: it is the excluded orphan described above.
