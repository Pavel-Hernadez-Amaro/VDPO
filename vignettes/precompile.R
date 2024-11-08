# move to the vignettes directory
old_dir <- setwd("vignettes")

# generate static Rmds
knitr::knit("source/_VDPO-01-introduction.Rmd", "VDPO-01-introduction.Rmd")
knitr::knit("source/_VDPO-02-vd-model.Rmd", "VDPO-02-vd-models.Rmd")

# restore
setwd(old_dir)
