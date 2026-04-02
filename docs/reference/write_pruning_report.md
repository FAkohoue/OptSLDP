# Write pruning statistics report

Write pruning statistics report

## Usage

``` r
write_pruning_report(pruning_stats, stats_file, summary_file = NULL)
```

## Arguments

- pruning_stats:

  A data.frame or data.table of step-by-step pruning counts as returned
  in the `pruning_stats` element of
  [`run_sldp()`](https://FAkohoue.github.io/OptSLDP/reference/run_sldp.md).

- stats_file:

  Path for the CSV output file.

- summary_file:

  Optional path for a plain-text summary file. If `NULL` (default) no
  summary is written.

## Value

Invisibly returns a named list with elements `stats_file` and
`summary_file`.
