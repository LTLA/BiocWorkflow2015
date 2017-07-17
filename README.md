# A Bioconductor workflow for DB analyses of ChIP-seq data

First we compile the R markdown file. 
We save the intermediates so that we can check the results against the expected values.

```r
# Execute in a R session
rmarkdown::render("workflow.Rmd", clean=FALSE)
```

We update the old result file; a simple `git diff` operation will reveal any differences between them.

```sh
cp workflow.knit.kmd workflow.md
```

Currently, alignment is turned off to ensure we get the same results upon re-analysis (as Rsubread changes quite regularly, and I don't won't to have to keep on checking it). 
However, we still need to check that the code is correct. 
To do that, we run it with alignment turned back on.

```sh
cd aligntest
cat ../markdown/workflow.Rmd | sed "s/remap <- FALSE/remap <- TRUE/" > workflow.Rmd
echo "rmarkdown::render('workflow.Rmd', clean=FALSE)" | R --no-save --vanilla
cd -
```

We also need to convert the R markdown into something that BioC can compile.
This is done using the 'toBioC.sh' script, which flips some flags and modifies the YAML metadata so that it is compatible with BioC's build system.

```
./toBioC.sh
```
