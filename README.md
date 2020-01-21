# Contents

1. [Getting started](#getting-started)
2. [Requirements](#requirements)
   1. [System libraries](#system-libraries)
   2. [R packages](#r-packages)
   3. [R packages installation guides](#r-packages-installation-guides)
      * [Installing R packages](#installing-r-packages)
      * [Building binary R packages from source](#building-binary-r-packages-from-source)
      * [Installing binary R packages](#installing-binary-r-packages)

# Getting started

```r
library(orynfect)
library(GloCR)
library(ggplot2)

WEATHER_DIR = "files/PHL"
RICE_SOS    = "files/CropCal/WORLD_PLANT_PK1_5.tif"
YEARS       = 2010:2018

infect <- function(disease=leafBlast, summary.fun=sum, field="severity", ...){
  infection <- disease(...)@d
  return(summary.fun(infection[,field]))
}

files.wth <- dir(WEATHER_DIR, pattern=".csv", full.names = TRUE)

rst.riceplant <- raster(RICE_SOS)
norice <- vector()

dat.infection <- vector()
for (i in 1:length(files.wth)){
  info.wth <- unlist(strsplit(basename(files.wth[i]),"_"))
  xy <- c(as.numeric(info.wth[3]),as.numeric(info.wth[4]))
  cell <- cellFromXY(rst.riceplant,xy)
  doy.emergence <- rst.riceplant[cell]
  
  if(is.na(doy.emergence)| cell %in% norice) {
    norice <- unique(c(norice,cell))
    next
  }
  dat.wth <- read.csv(files.wth[i], stringsAsFactors = FALSE)
  dat.wth$date <- as.Date(dat.wth$date)
  colnames(dat.wth) <- c("date", "tmax", "tmin", "prec", "srad", "rhmax","rhmin", "wind.morningmax", "wind.daymax", "wind.avg")
  dat.wth$tavg <- (dat.wth$tmax + dat.wth$tmin)/2

  wth <- dat.wth

  for(yr in YEARS){
    message("Cell-",cell, "_Year-",yr)
    infection <- data.frame(cell, year=yr, season="main", audpc=infect(wth=wth, crop.estabdate=dateFromDoy(doy.emergence,yr)))
    dat.infection <- rbind(dat.infection, infection)
  }
  
}

dat.infection <- cbind(dat.infection, xyFromCell(rst.riceplant, dat.infection$cell)) 
dat.infection$audpc.class <- cut(dat.infection$audpc,breaks=seq(0,450,length.out = 9))
ggplot()+ geom_tile(data = dat.infection, aes(x=x,y=y, fill=audpc.class )) +facet_wrap(~year) + scale_fill_brewer(palette = "Reds") + theme_dark() + theme(legend.position = "bottom")
ggsave("files/sample.leafblast.png")
```

# Requirements

## System libraries

Whether you are on Windows, Mac, or Linux you will need to install the required system libraries:

1. [GDAL](https://gdal.org/)

## R packages

Below are the R packages needed to run the R script:

1. `orynfect` - **build**-and-**install** from source
2. `manipulateR`
3. `RODBC`
4. `bitops`
5. `RCurl`
6. `GloCR`
7. `raster`
8. `rgdal`

## R package installation guides

### Installing R packages

1. In the terminal, issue the statement:
   ```bash
   R
   ```

   It will output something like:
   ```
   R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night"
   Copyright (C) 2019 The R Foundation for Statistical Computing
   Platform: x86_64-apple-darwin15.6.0 (64-bit)

   R is free software and comes with ABSOLUTELY NO WARRANTY.
   You are welcome to redistribute it under certain conditions.
   Type 'license()' or 'licence()' for distribution details.

   Natural language support but running in an English locale

   R is a collaborative project with many contributors.
   Type 'contributors()' for more information and
   'citation()' on how to cite R or R packages in publications.

   Type 'demo()' for some demos, 'help()' for on-line help, or
   'help.start()' for an HTML browser interface to help.
   Type 'q()' to quit R.

   >
   ```
2. Issue the statement:
   ```
   install.packages("manipulateR")
   ```

   Replace `manipulateR` with the package name you want to install.


### Building binary R packages from source

1. Download the latest copy of the R raw package - [georyza-master.zip](https://github.com/jaunario/georyza/archive/master.zip)
2. Unzip the `master.zip` file
3. In the terminal, go to the unzip directory. Issue the statement:
   ```bash
   cd /path/to/georyza
   ```
4. Issue the statement:
   ```bash
   R CMD build orynfect
   ```
   The command statement will output a binary tarball file like so: `orynfect_x.x.x.tar.gz`

### Installing binary R packages

1. Download a copy of [orynfect_0.0.1.tar.gz](#) and [orysat_0.0.3.tar.gz](#) packages. Alternatively, you can [build your own binary packages](#building-binary-packages).

2. In the terminal, go to the download directory. Issue the statement:
   ```bash
   cd /path/to/the/download/dir
   ```

3. Issue the statement:
   ```bash
   R CMD INSTALL orynfect_0.0.1.tar.gz
   ```



### Welcome to GitHub Pages

You can use the [editor on GitHub](https://github.com/jaunario/georyza/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/jaunario/georyza/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
