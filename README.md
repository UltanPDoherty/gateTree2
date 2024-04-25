gateTree
================
Ultán P. Doherty
2024-04-25

``` r
library(healthyFlowData)
library(ggplot2)
library(devtools)
load_all()
```

``` r
data(hd)
hfd1 <- hd.flowSet[[1]]@exprs

GGally::ggpairs(hfd1, progress = FALSE)
```

![](README_files/figure-gfm/hfd1_setup-1.png)<!-- -->

``` r
plusminus <- rbind(
  "CD4+_T" = c(+1, -1, +1, -1),
  "CD8+_T" = c(-1, +1, +1, -1),
  "B"      = c(-1, -1, -1, +1)
)
colnames(plusminus) <- colnames(hfd1)
plusminus
```

    ##        CD4 CD8 CD3 CD19
    ## CD4+_T   1  -1   1   -1
    ## CD8+_T  -1   1   1   -1
    ## B       -1  -1  -1    1

``` r
hfd1_gatetree <- gatetree(hfd1, plusminus, min_scaled_bic_diff = 50)
```

``` r
hfd1_gatetree$tree_plot + 
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_continuous(expand = c(0.1, 0.1))
```

![](README_files/figure-gfm/tree_plot-1.png)<!-- -->

``` r
GGally::ggpairs(hfd1, progress = FALSE,
                aes(colour = as.factor(1 + hfd1_gatetree$labels))) +
  ggokabeito::scale_colour_okabe_ito(order = c(9, 1, 2, 3)) +
  ggokabeito::scale_fill_okabe_ito(order = c(9, 1, 2, 3))
```

![](README_files/figure-gfm/ggpairs-1.png)<!-- -->
