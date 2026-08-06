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
        BENCHMARK/COMPATIBILITY
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

Measured on 2026-07-26 against the UCSC `linux.x86_64` binaries and
`NetFilterNonNested.perl` from `hillerlab/GenomeAlignmentTools`.

## Environment

- CPU: AMD Ryzen 7 5700X, 8 cores / 16 threads
- OS: Linux 6.14.5-100.fc40.x86_64
- Compiler: rustc 1.94.1, `cargo build --release`
- Reference binaries: `chainCleaner`, `chainNet`, `chainSort` from
  <https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/>, plus
  `NetFilterNonNested.perl`
- Timer and memory: GNU `time -v`, one run per row

## Inputs

| Input | Detail |
|:--|:--|
| Chains | `hg38.HLrhiRex1.filled.all.chain.gz`, score-sorted with `chaintools sort` (10.3 s, excluded from all timings) |
| Score-sorted chain | 722,463,293 B, 3,049,608 chains, 608 reference sequences, 91 query sequences, 2 metadata lines, SHA-256 `ad48641c…` |
| Reference | `hg38.p14.2bit` from UCSC `goldenPath/hg38/bigZips/p14`, 711 sequences |
| Query | `GCA_041825415.1` (*Rhinolophus rex*, ASM4182541v1) converted with `faToTwoBit`, sequence names reduced to the accession without version, 149 sequences |
| Sizes | `twoBitInfo` output for both 2bit files, given to both implementations |
| Options | `-linearGap=loose` and otherwise stock defaults |

The downloaded chain file is not score-sorted, so both `chainNet` and the Rust
port reject it; sorting once up front is what makes the two pipelines
comparable. The patched hg38 is required because 255 of the 608 reference
sequences in the chain file are `_fix`/`_alt` scaffolds absent from the plain
`hg38.2bit`.

## Whole-genome results

### NET supplied to both (`-net` / `--net`)

The NET is `chainNet -minScore=0` followed by
`NetFilterNonNested.perl -minScore1 3000`, produced once and handed to both.

| Command | Wall | Peak RSS | Speedup |
|:--|--:|--:|--:|
| Kent `chainCleaner` | 55.97 s | 10.26 GiB | 1.00× |
| `chainc` | 24.61 s | 10.96 GiB | 2.27× |

For context, producing that NET costs `chainNet` 18.74 s / 3.63 GiB plus
`NetFilterNonNested.perl` 24.36 s / 0.46 GiB.

### No NET given: each implementation nets the chains itself

Kent shells out to `chainNet` and the Perl filter through a temporary file;
`chainc` builds the same NET in memory through `netools`, reference side
only, with the non-nested filter fused into construction.

| Command | Wall | Peak RSS | Speedup |
|:--|--:|--:|--:|
| Kent `chainCleaner` (nets internally) | 99.16 s | 10.26 GiB | 1.00× |
| `chainc` (nets in memory) | 28.31 s | 10.85 GiB | 3.50× |

### Thread counts

| Threads | Wall | Peak RSS |
|--:|--:|--:|
| 1 | 30.36 s | 10.22 GiB |
| 16 | 27.85 s | 10.76 GiB |

Both runs are byte-identical to each other and to the run above. Sequence
loading and chain parsing dominate, so extra threads buy little here.

## Output parity

Both implementations emit 3,050,726 chains: the 3,049,608 input chains plus
1,118 new chains for removed suspects.

| Output | Result |
|:--|:--|
| `removed.bed` | byte-identical, 1,118 records, 195,765 B |
| `newChainIDDict` | byte-identical, 1,118 records |
| Chain metadata | byte-identical, both `##matrix` and `##gapPenalties` lines |
| Chain records | identical multiset, 3,050,726 records, 722,560,346 B in both |
| Chain scores | identical score at every output position |
| Generated NET | byte-identical to `chainNet` + Perl: 113,168,312 B, MD5 `4bdbf67807ff057cfa2f9c9f96a1c977` |

Each implementation also produces exactly the same output with and without
`--net`, so the netting path changes nothing downstream.

The chain files are not byte-identical, and cannot be: `chainSort` compares only
the score and sorts through C `qsort`, which is not stable. 98.82% of records
fall in equal-score groups (96,817 groups, largest 1,629 records), and every
position where the two files disagree lies inside such a group. Every record
with a unique score appears at the same position with the same chain ID in both
files.

## Smaller runs

`chr21` + `chr22` subset, 67,929 input chains, both modes:

| Command | NET supplied | Nets internally |
|:--|--:|--:|
| Kent `chainCleaner` | 5.78 s | 6.99 s |
| `chainc` | 5.79 s | 5.93 s |

Outputs agree as above (44 removals, BED and dictionary byte-identical, NET
byte-identical). At this size the run is dominated by loading 43 query
chromosomes and 2 reference chromosomes — 3.5 GiB in both implementations — so
there is no speedup to be had.

`bench/diff_fixture.sh` compares the two implementations on the unit fixtures,
where output including chain records is byte-identical for a supplied NET, for
`-doPairs`, and for the netting path.

## Reproducing

```bash
bench/diff_fixture.sh single bench/fix/reference.2bit bench/fix/query.2bit \
  -net=tests/fixtures/single.net -linearGap=loose

bench/diff_real.sh full-net bench/sorted.chain bench/hg38.p14.2bit \
  bench/rhiRex1.2bit bench/t14.sizes bench/q.sizes net
bench/diff_real.sh full-nonet bench/sorted.chain bench/hg38.p14.2bit \
  bench/rhiRex1.2bit bench/t14.sizes bench/q.sizes nonet
```

`bench/` holds roughly 11 GiB of downloaded and intermediate data, including
`rhiRex1.fa` (2.1 GiB) and the unsorted `input.chain`, both only needed to
rebuild `rhiRex1.2bit` and `sorted.chain`.
