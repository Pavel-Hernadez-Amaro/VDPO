# move to the vignettes directory
old_dir <- setwd("vignettes")

# generate static Rmds
knitr::knit("source/_VDPO-01-introduction.Rmd", "VDPO-01-introduction.Rmd")
knitr::knit("source/_VDPO-02-vd-model.Rmd", "VDPO-02-vd-models.Rmd")
knitr::knit("source/_VDPO-03-pofd.Rmd", "VDPO-03-pofd.Rmd")
knitr::knit("source/_VDPO-04-vd-mfpca.Rmd", "VDPO-04-vd-mfpca.Rmd")

# restore
setwd(old_dir)
