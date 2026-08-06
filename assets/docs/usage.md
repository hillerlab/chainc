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
        USER GUIDE
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


```text
chainc \
  --chains input.chain \
  --net input.net \
  --reference reference.2bit \
  --query query.2bit \
  --output output.chain \
  --removed-bed removed.bed \
  --linear-gap loose
```

Sequence inputs may be `.2bit`, FASTA, or gzip-compressed FASTA. Run
`chainc --help` for all thresholds and auxiliary outputs.

## Install

```bash
cargo build --release
```

The binary is `target/release/chainc`.

## Outputs

| File | Content |
|:--|:--|
| `--output` | All chains, sorted by descending score, with original metadata lines preserved |
| `--removed-bed` | BED9 record per removed suspect: breaking/broken chain IDs and header scores, suspect local score, left/right ratios; RGB `0,0,153` for a single suspect, `0,100,255` for a pair |
| `--new-chain-id-dict` | One `NEW_CHAIN_ID<TAB>ORIGINAL_BREAKING_CHAIN_ID` line per removal |
| `--suspect-data-file` | Dry-run table of every single suspect; no removal happens and pair processing is disabled |

## Netting without an input NET

`--net` is optional. When it is omitted the chains are netted in memory, which
reproduces what the original program does by shelling out to

```text
chainNet -minScore=0 in.chain t.sizes q.sizes tmp.net.raw /dev/null
NetFilterNonNested.perl tmp.net.raw -minScore1 3000 > tmp.net
```

`netools` 0.0.2 ports both steps and fuses them, so only the reference side is
built, the raw NET is never materialized, and no temporary file or subprocess is
involved:

```text
chainc \
  --chains input.chain \
  --reference reference.2bit \
  --query query.2bit \
  --output output.chain \
  --removed-bed removed.bed \
  --linear-gap loose
```

Generated NET sections follow the order of the sequence sizes, exactly as
`chainNet` follows the order of its `tSizes` file. Sizes come from the sequence
files themselves — from the index only for `.2bit` — or from
`--reference-sizes` / `--query-sizes` when a specific order or subset is wanted.
As with `chainNet`, the input chains must be sorted by descending score.
`--dump-net` writes the NET that was used.

## Options

| Long | Short | Default | Meaning |
|:--|:-:|:--|:--|
| `--chains` | `-c` | required | Input chain file (`.chain` or `.gz`) |
| `--net` | `-n` | none | Hierarchical NET; when omitted, chains are netted in memory |
| `--reference-sizes` | `-R` | sequence order | Ordered reference chromosome sizes for the generated NET |
| `--query-sizes` | `-Q` | sequence order | Ordered query chromosome sizes for the generated NET |
| `--dump-net` | `-N` | none | Write the NET used for discovery (hidden) |
| `--reference` | `-r` | required | Reference genome, `.2bit`, FASTA, or gzipped FASTA |
| `--query` | `-q` | required | Query genome, `.2bit`, FASTA, or gzipped FASTA |
| `--output` | `-o` | required | Sorted output chain file |
| `--removed-bed` | `-b` | required | BED9 output of removed suspects |
| `--linear-gap` | `-g` | required | `loose`, `medium`, or a Kent-style gap file |
| `--score-scheme` | `-s` | AXT DNA matrix | Blastz/LASTZ-format substitution matrix |
| `--fold-threshold` | `-F` | `0.0` | Minimum `fill_global / suspect_local` ratio |
| `--lr-fold-threshold` | `-L` | `2.5` | Minimum left/right fill ratio for a single suspect |
| `--lr-fold-threshold-pairs` | `-P` | `10.0` | Minimum left/right fill ratio for a pair |
| `--max-suspect-bases` | `-B` | `i32::MAX` | Maximum aligned bases of an accepted suspect |
| `--max-suspect-score` | `-M` | `100000` | Maximum local score of an accepted suspect |
| `--min-broken-chain-score` | `-S` | `50000` | Minimum broken-chain header score |
| `--min-lr-gap-size` | `-G` | `0` | Minimum left and right enclosing gap size |
| `--do-pairs` | `-p` | off | Also test adjacent pairs of suspects |
| `--max-pair-distance` | `-D` | `10000` | Maximum reference distance between paired suspects |
| `--new-chain-id-dict` | `-i` | none | Write `NEW_CHAIN_ID<TAB>ORIGINAL_BREAKING_CHAIN_ID` |
| `--suspect-data-file` | `-d` | none | Dry-run table of single suspects; disables removal and pairs |
| `--threads` | `-t` | all CPUs | Worker threads for parallel parsers |
| `--strict` | `-e` | off | Reject ranking anomalies and zero-local-score candidates |
| `--dump-candidate-order` | `-x` | none | Write the compatibility-ordered candidate trace (hidden) |

A suspect is removed only when every condition passes: both flanking fills
score at least `--lr-fold-threshold` (or `--lr-fold-threshold-pairs` for pairs)
against the suspect's local score, the complete fill scores at least
`--fold-threshold`, the suspect's local score and aligned bases stay under
`--max-suspect-score` / `--max-suspect-bases`, the broken chain's header score
is at least `--min-broken-chain-score`, and both enclosing gaps are at least
`--min-lr-gap-size` bases. Pairs must additionally share depth, broken chain,
and middle gap, and be at most `--max-pair-distance` apart.

## Compatibility

The default behavior follows the original program's scoring and mutation
semantics:

- Chain header scores are not updated during removal, so a chain that breaks one
  chain and is later broken by another can be tested against a stale score.
- Scoring a subset that spans a whole chain rewrites that chain's score, which is
  what `chainSubsetOnT` plus `getChainScore` do when the requested interval
  covers the chain. Later thresholds, BED records, and the written chain all see
  the new value.
- Break candidates, breaking-chain groups, new chain IDs, and BED order follow
  Kent's hash layout, reproduced from `hashString`, the 4096-bucket initial
  table, its expansion rule, and `hashElListHash`.

Output is deterministic. Ties are ordered by output ordinal and then chain ID
because the C `qsort` inside `chainSort` does not define equal-score order.

`bench/diff_fixture.sh` checks chain, BED, and ID-dictionary output byte for byte
against the UCSC binary, both for a supplied NET and for the netting path.
