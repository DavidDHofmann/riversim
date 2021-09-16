
# riversim <img src="man/figures/riversim.png" align="right" width="150" height="150"/>

`riversim` is an R-package that allows you to simulate river networks.
The algorithm was originally developed by a [reddit
user](https://www.reddit.com/r/proceduralgeneration/comments/ftgbgo/ive_been_working_on_an_algorithm_to_efficiently/).

## Installation

To install `riversim` you need to have the package `devtools` installed.
You can simply do:

``` r
library(devtools)
install_github("DavidDHofmann/riversim")
```

## Example

Here is a simple example of a river network simulation

``` r
library(riversim)

# Generate the river network
river <- rivernetwork(1000)

# Visualize it
plot(river, col = "cornflowerblue", border = NA)
axis(1)
axis(2)

# The network can also be returned as a raster
river <- rivernetwork(1000, raster = T)
plot(river, col = c("white", "cornflowerblue"))
```
