[![Contributors][contributors-shield]][contributors-url]
[![MPL2.0 License][license-shield]][license-url]

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <h3 align="center">TidyComb</h3>

  <p align="center">
    Tidy drug combination high-throughput screening data for DrugComb portal
    <br />
    <a href="https://irscope.shinyapps.io/chloroplot/"><strong>Chloroplot</strong></a>
    <p> An R package for parsing and visualizing organelle genomes</p>
    <br />
    <a href="https://github.com/shuyuzheng/Chloroplot/issues">Report Bug</a>
    ·
    <a href="https://github.com/shuyuzheng/Chloroplot/issues">Request Feature</a>
  </p>
</p>


<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Package](#about-the-package)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
* [Reference](#reference)
* [License](#license)
* [Contact](#contact)

## About the Package 

The R package **Chloroplot** wrapped the functions for visualizing the organelle genomes.
We also provide a web application at [https://irscope.shinyapps.io/chloroplot/](https://irscope.shinyapps.io/chloroplot/) for users who prefere interactive GUI.

### Built with

* [R](https://www.r-project.org/) 3.6
* [RStudio](https://www.rstudio.com/) 1.2

## Getting Start

### Prerequistes

**Chloroplot** is a package written with [R](https://www.r-project.org/) programming language. Please make sure that you have installed [R project](https://www.r-project.org/) ( >= V3.6) on your local machine. You can download and install R from [CRAN](https://cran.r-project.org/mirrors.html).

### Installation

You can install **Chloroplot** from this [GitHub repository](https://github.com/shuyuzheng/Chloroplot):

1. Install and load the [devtools](https://github.com/hadley/devtools) package. You can do this from [CRAN](https://cran.r-project.org/). Invoke R and then type

```
install.packages("devtools")
library(devtools)
```

2. To install **Chloroplot** from [GitHub](https://github.com/), you'd type:

```
devtools::install_github("DrugComb/Chloroplot")
```

## Usage

There are 3 main functions in Chloroplot for visualizing the organelle genomes: "PlotTab", " PlotMitGenome"and "PlotPlastidGenome". Following are two examples for utality:

```
library(chloroplot)
# Plot Ppchloroplast genome.

# 1. Parsing the GenBank file.
#"EU549769" is the GenBank accession for Guizotia abyssinica chloroplast 
t <- PlotTab(gbfile = "EU549769")
# 2. Generate plot
PlotPlastidGenome(t) # The plot will be saved in a pdf file under your work directory.
PlotPlastidGenome(t, save = FALSE) # The plot will be shown in the "plot" panel if you are using Rstudio.

# Plot mitochondrion genome.

# 1. Parsing the GenBank file
# "NC_012920.1" is the GenBank accession for Homo sapiens mitochondrion
t <- PlotTab(gbfile = "NC_012920.1")
# 2. Generate plot
PlotMitGenome(t) # The plot will be saved in a pdf file under your work directory.
PlotMitGenome(t, save = FALSE) # The plot will be shown in the "plot" panel if you are using Rstudio.
```

For further details about manipulating the visualizations. Please check the documentations for the functions, by typing following commands in the R console.

```
help(PlotMitGenome) # or ?PlotMitGenome
help(PlotPlastidGenome) # or ?PlotPlastidGenome
```

## Reference

[1]: Zheng, S.; Poczai, P.; Hyvönen, J.; Tang, J.; Amiryousefi, A. Chloroplot: An Online Program for the Versatile Plotting of Organelle Genomes. Front. Genet. 2020, 11. https://doi.org/10.3389/fgene.2020.576124.

## License

Distributed under the Mozilla Public License 2.0

## Contact

Shuyu Zheng - shuyu.zheng@helsinki.fi
ali.amiryousefi@helsinki.fi - Ali Amiryousefi

Project Link: https://github.com/shuyuzheng/Chloroplot

## Acknowledgements
* [Img Shields](https://shields.io)
* [Choose an Open Source License](https://choosealicense.com)

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/badge/contributors-1-orange.svg?style=flat-square
[contributors-url]: https://github.com/shuyuzheng/Chloroplot/graphs/contributors
[license-shield]: https://img.shields.io/badge/license-MPL--2.0-blue.svg
[license-url]: https://choosealicense.com/licenses/mpl-2.0
