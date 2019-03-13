Class07 Functions Part 2
================
Gidsela Luna
January 30, 2019

Functions Revisted
==================

Load(i.e **source**) our rescale() function from last time

``` r
source("http://tinyurl.com/rescale-R")
```

``` r
rescale(1:5)
```

    ## [1] 0.00 0.25 0.50 0.75 1.00

``` r
#rescale(c(1:5,"string"))

#Error in x - rng[1] : non-numeric argument to binary operator
```

### We want to make this function more robust to these types of errors

``` r
#rescale2(c(1:5,"string"))

#Error: Input x should be numeric
```

``` r
is.numeric(1:5)
```

    ## [1] TRUE

``` r
is.numeric("string")
```

    ## [1] FALSE

``` r
is.numeric(c(1:5,"string"))
```

    ## [1] FALSE

``` r
!is.numeric(c(1:5,"string"))
```

    ## [1] TRUE

``` r
!is.numeric(1:5)
```

    ## [1] FALSE

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)==TRUE
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
which(is.na(x) & is.na(y) == TRUE)
```

    ## [1] 3

``` r
sum(is.na(x) & is.na(y)==TRUE)
```

    ## [1] 1

``` r
both_na<- function(x,y){
  #check for the number of NA elements in the vector 
  sum(is.na(x) & is.na(y)==TRUE)
}
```

``` r
both_na(x,y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
y3 <- c(1,NA,NA,NA,NA)

#what will this return
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y) == TRUE: longer object length is not a
    ## multiple of shorter object length

    ## [1] 3

``` r
# it is recycling
both_na(x,y3)
```

    ## Warning in is.na(x) & is.na(y) == TRUE: longer object length is not a
    ## multiple of shorter object length

    ## [1] 4

``` r
#both_na2(x,y2)

#Error: Input x and y should be vectors of the same length
```

ˆ

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
both_na3(x,y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3
