# Grid adder for dataframes

It prepared the partially observed data to be inputed in the `ffpo`
function. This function should only be used when the
`bidimensional_grid` parameter of the `ffpo` function is `FALSE`.

## Usage

``` r
add_grid(df, grid)
```

## Arguments

- df:

  `data.frame` object to which the grid will be added.

- grid:

  Grid vector.

## Value

`data.frame` with the grid added.

## See also

[`ffpo`](https://pavel-hernadez-amaro.github.io/VDPO/reference/ffpo.md)
