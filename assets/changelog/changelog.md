<p align="center">
  <p align="center">
    <img width=100 align="center" src="../figures/chainc.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="../figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="../figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="../figures/hillerlab-light.png"
    >
  </picture>
</p>

  <span>
    <h1 align="center">
        chainc
    </h1>
  </span>

  <span>
    <h1 align="center">
        CHANGELOG
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/chainc" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/chainc?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://genome.ucsc.edu/goldenpath/help/chain.html">chains</a> .
        <a href="https://github.com/alejandrogzi/chainc/blob/master/assets/docs/usage.md">usage</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>


## 0.0.3 - 2026-11-20

### Changed

- Bumped `chaintools` from 0.0.12 to 0.0.14 and `netools` from 0.0.3 to 0.0.4.
  Sequences are now loaded with the one-pass 2bit loader of `chaintools`
  0.0.14, which reads version 1 `.2bit` files (64-bit offsets, written by
  `faToTwoBit -long` for genomes above 4 Gb). Earlier releases stopped on them
  with `Unsupported version: Versions larger than 0 are not supported`.
- No source change. Output is byte-identical to 0.0.1 and 0.0.2 (chain,
  `removed.bed`, `newChainIDDict`) on version 0 `.2bit` and FASTA input.

### Known limitation

- Without `--reference-sizes` / `--query-sizes` and without `--net`, sizes are
  read from the `.2bit` index through `netools`, which still rejects version 1
  files. Pass the size files (or a NET) for such genomes.

## 0.0.1 - 2026-08-06

First release: a Rust port of UCSC `chainCleaner`.

### Added

- Kent-compatible candidate discovery from a hierarchical NET, including the
  span index, lower-chain-ID overlap filtering, and pair construction
- Exact `chainSubsetOnT` behavior, internal block removal, and iterative
  neighbor-boundary updates within each breaking-chain group
- Kent-compatible global and local scoring over `.2bit` / FASTA / gzipped-FASTA
  sequences, with per-revision metrics caching
- Break-candidate and removal order reproduced from Kent's hash layout
  (`hashString`, 4096-bucket table, expansion, `hashElListHash`)
- In-memory netting through `netools` when `--net` is omitted, fusing
  `chainNet -minScore=0` and `NetFilterNonNested.perl -minScore1 3000` into one
  reference-side pass with no temporary files or subprocesses
- CLI for every threshold, pair mode, dry-run `--suspect-data-file`, strict
  mode, `--new-chain-id-dict`, and hidden `--dump-net` / `--dump-candidate-order`
- Deterministic output at any thread count, with a documented tie-breaker for
  equal scores
- Differential fixtures in `bench/diff_fixture.sh`

### Performance

Byte-identical `removed.bed` (1,118 records), `newChainIDDict`, chain metadata,
and generated NET against the UCSC binaries on a whole-genome run; chain records
identical as a multiset with identical scores everywhere. Wall-clock speedups:
2.27x with a supplied NET, 3.50x when netting in memory. Full numbers in
[bench.md](../docs/bench.md).
