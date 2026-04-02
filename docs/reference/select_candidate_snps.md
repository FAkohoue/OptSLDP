# Select candidate SNPs by screening threshold

Applies one of three selection modes to a screening statistics table.

## Usage

``` r
select_candidate_snps(
  stats_dt,
  mode = c("A", "B", "C"),
  pval_threshold = NULL,
  z_threshold = NULL,
  pve_threshold = NULL,
  logic = c("AND", "OR")
)
```

## Arguments

- stats_dt:

  `data.table` from
  [`compute_screening_stats()`](https://FAkohoue.github.io/OptSLDP/reference/compute_screening_stats.md)
  for a **single trait**.

- mode:

  `"A"` (default), `"B"`, or `"C"`.

- pval_threshold:

  P-value upper bound (modes A / C).

- z_threshold:

  \|z-score\| lower bound (modes B / C).

- pve_threshold:

  PVE lower bound (modes B / C).

- logic:

  `"AND"` (default) or `"OR"`.

## Value

Character vector of selected SNP IDs.

## Details

- Mode A:

  P-value only. Requires `pval_threshold`.

- Mode B:

  Effect-size criteria. Requires `z_threshold` and/or `pve_threshold`.

- Mode C:

  Combination of A and B criteria.
