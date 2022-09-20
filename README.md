-   [Bayesian-Consensus-Clustering](#bayesian-consensus-clustering)
    -   [Setup](#setup)
    -   [Result from summary folder](#result-from-summary-folder)
    -   [Results from graphic output
        folder](#results-from-graphic-output-folder)
        -   [Trajectories for each
            marker](#trajectories-for-each-marker)

# Bayesian-Consensus-Clustering

## Setup

There are two main folders The summary folder contains code that gives
model summary only. The graphic output folder visualizes the results. To
run the code, first navigate to either `summary` or `graphich output`
folder and set current folder as working directory.

## Result from summary folder

Results from exmaple

``` r
> fit.BCC$summary.stat$PPI
                  [,1]        [,2]
mean        0.79474912  0.20525088
sd          0.04368524  0.04368524
2.5%        0.70057013  0.13256579
97.5%       0.86743421  0.29942987
geweke.stat 0.92788809 -0.92788809

> fit.BCC$summary.stat$ALPHA
                   [,1]       [,2]        [,3]
mean         0.80260557 0.50588271  0.87723123
sd           0.04409383 0.00546498  0.05000172
2.5%         0.73086562 0.50011619  0.78506714
97.5%        0.87966656 0.51894468  0.98112848
geweke.stat -0.44096846 0.28068139 -0.24425527

> fit.BCC$summary.stat$GA
[[1]]
, , 1

                   [,1]      [,2]
mean        -0.40601148 1.5927662
sd           0.06745642 0.1047098
2.5%        -0.53719087 1.4105057
97.5%       -0.26418628 1.8241435
geweke.stat -1.77950969 1.5074619

, , 2

                    [,1]        [,2]
mean         0.004370976 0.014538832
sd           0.002038668 0.003148562
2.5%         0.000826820 0.008778884
97.5%        0.008498853 0.020473345
geweke.stat -3.021475187 0.775005515


[[2]]
, , 1

                  [,1]        [,2]
mean        4.70529869  5.59804230
sd          0.05915515  0.03607394
2.5%        4.58360223  5.53053336
97.5%       4.80947449  5.66298393
geweke.stat 0.22801624 -0.75583503

, , 2

                    [,1]         [,2]
mean        -0.011684530 -0.002390531
sd           0.001637444  0.001201627
2.5%        -0.014946694 -0.004739064
97.5%       -0.008976688 -0.000243194
geweke.stat  0.470934639 -0.932315679


[[3]]
, , 1

                  [,1]      [,2]
mean        -1.5606782 1.7970701
sd           0.4618175 0.6398377
2.5%        -2.4133218 0.7165736
97.5%       -0.5889110 3.3114258
geweke.stat  0.8325084 1.1956285

, , 2

                   [,1]        [,2]
mean        0.030584706  0.01236938
sd          0.015885539  0.01737053
2.5%        0.005538517 -0.02911865
97.5%       0.070120505  0.04451236
geweke.stat 0.274449168  0.67597695


> fit.BCC$summary.stat$SIGMA.SQ.U
[[1]]
, , 1

                     [,1]
mean         3.727492e-05
sd           8.293761e-06
2.5%         2.584153e-05
97.5%        5.369400e-05
geweke.stat -2.112115e-01

, , 2

                    [,1]
mean        5.217424e-05
sd          9.666406e-06
2.5%        3.673481e-05
97.5%       7.113944e-05
geweke.stat 6.282917e-01


[[2]]
, , 1

                     [,1]
mean         3.409460e-05
sd           6.500162e-06
2.5%         2.403157e-05
97.5%        4.874368e-05
geweke.stat -2.053133e+00

, , 2

                     [,1]
mean         2.319546e-05
sd           3.611774e-06
2.5%         1.643193e-05
97.5%        3.121602e-05
geweke.stat -1.554764e-01


[[3]]
, , 1

                    [,1]
mean        4.102649e-05
sd          8.384566e-06
2.5%        2.842640e-05
97.5%       5.941488e-05
geweke.stat 4.752748e-01

, , 2

                    [,1]
mean        6.359277e-05
sd          1.745721e-05
2.5%        3.921035e-05
97.5%       1.040450e-04
geweke.stat 1.242202e+00


> fit.BCC$summary.stat$SIGMA.SQ.E
[[1]]
                  [,1]       [,2]
mean        0.37393750 0.37393750
sd          0.01775766 0.01775766
2.5%        0.34146760 0.34146760
97.5%       0.40594660 0.40594660
geweke.stat 1.96269393 1.96269393

[[2]]
                   [,1]        [,2]
mean        0.084539874 0.084539874
sd          0.003916127 0.003916127
2.5%        0.078067818 0.078067818
97.5%       0.093544052 0.093544052
geweke.stat 0.581838511 0.581838511

[[3]]
NULL

> fit.BCC$alpha.adjust
[1] 0.4571463
```

## Results from graphic output folder

### Trajectories for each marker

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)
