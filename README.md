# CloneDeMix

This is an R package designed for inferring subclonal structure of tumors on copy number aberrations.

To install this package, please copy and paste following codes into your R session:

1. install.packages("devtools")
2. library(devtools)
3. install_github("AshTai/CloneDeMix")

## Example
```{r settings, include = FALSE}
library(CloneDeMix)
data("ESCC_chr1")
res <- CloneDeMix(tumor=ESCC_chr1$tumor, normal=ESCC_chr1$normal,threshold = 10^-5, iterC = 10^3,
  CNVstate = c(0:10), method = "aic")
head(res$CNV); head(res$MCP)
```
