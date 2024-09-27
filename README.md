[![Rust](https://github.com/stuart-lab/sparx/actions/workflows/rust.yml/badge.svg)](https://github.com/stuart-lab/sparx/actions/workflows/rust.yml)

# sparx: disk-based sparse matrix tools

`sparx` is a memory-efficient tool for working with sparse matrices in the Matrix Market format on-disk.

## Usage

Subset a matrix on-disk:

```
sparx subset -i <matrix.mtx> --rows <row_index.txt> --cols <col_index.txt> -o <matrix_subset.mtx>
```

Compute matrix statistics (nonzero count, sum, mean, variance, standard deviation, min, max) for each row and column:

```
sparx stats -i <matrix.mtx> -o <filename>
```

## Installation

To compile, clone the git repo and run `cargo install`:

```
git clone git@github.com:stuart-lab/sparx.git
cd sparx; cargo install --path .
```