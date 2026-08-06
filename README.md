<p align="center">
  <p align="center">
    <img width=100 align="center" src="./assets/figures/chainc.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="./assets/figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="./assets/figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="./assets/figures/hillerlab-light.png"
    >
  </picture>
</p>

  <span>
    <h1 align="center">
        chainc
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/chainc" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/chainc?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> Remove chain-breaking alignments using chain/net files </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://genome.ucsc.edu/goldenpath/help/chain.html">chains</a> .
        <a href="https://github.com/alejandrogzi/chainc/blob/master/assets/docs/usage.md">usage</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

> [!IMPORTANT]
> `chainc` is a Rust implementation of UCSC `chainCleaner`. It discovers
> chain-breaking alignments from a hierarchical net, evaluates them with
> Kent-compatible DNA and linear-gap scores, removes accepted internal blocks,
> and writes sorted chains plus BED diagnostics.

---

# Benchmarks

See [bench.md](assets/docs/bench.md) for full details.

### .net supplied to both (`-net` / `--net`)

The NET is `chainNet -minScore=0` followed by
`NetFilterNonNested.perl -minScore1 3000`, produced once and handed to both.

<div align="center">

| Command | Wall | Peak RSS | Speedup |
|:--|--:|--:|--:|
| Kent `chainCleaner` | 55.97 s | 10.26 GiB | 1.00× |
| `chainc` | 24.61 s | 10.96 GiB | 2.27× |

</div>

For context, producing that NET costs `chainNet` 18.74 s / 3.63 GiB plus
`NetFilterNonNested.perl` 24.36 s / 0.46 GiB.

### No .net given: each implementation nets the chains itself

Kent shells out to `chainNet` and the Perl filter through a temporary file;
`chaincleaner` builds the same NET in memory through `netools`, reference side
only, with the non-nested filter fused into construction.

<div align="center">

| Command | Wall | Peak RSS | Speedup |
|:--|--:|--:|--:|
| Kent `chainCleaner` (nets internally) | 99.16 s | 10.26 GiB | 1.00× |
| `chainc` (nets in memory) | 28.31 s | 10.85 GiB | 3.50× |

</div>

### Output parity

Both implementations emit 3,050,726 chains: the 3,049,608 input chains plus
1,118 new chains for removed suspects.

<div align="center">

| Output | Result |
|:--|:--|
| `removed.bed` | byte-identical, 1,118 records, 195,765 B |
| `newChainIDDict` | byte-identical, 1,118 records |
| Chain metadata | byte-identical, both `##matrix` and `##gapPenalties` lines |
| Chain records | identical multiset, 3,050,726 records, 722,560,346 B in both |
| Chain scores | identical score at every output position |
| Generated NET | byte-identical to `chainNet` + Perl: 113,168,312 B, MD5 `4bdbf67807ff057cfa2f9c9f96a1c977` |

</div>
