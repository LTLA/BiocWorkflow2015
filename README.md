# A Bioconductor workflow for DB analyses of ChIP-seq data

Compilation of the workflow is quite involved, so an explanation is provided here.
The first step is to convert the figure tags into figure numbers. 
This is because markdown doesn't have particularly good support for figure labelling - at least, not in a flexible manner that works for both LaTeX and HTML. 
Instead, we hard-code the figure numbers into the markdown text, so that it makes sense in all downstream applications. 
Figure numbers are also hardcoded into the captions, which means that you need to specify empty captions within LaTeX.

```
cd markdown
./refigure.sh
cd -
```

The second step is to compile the R markdown file. 

```
cd markdown
echo "knitr::knit('workflow.Rmd')" | R --no-save --vanilla
cd -
```

Currently, alignment is turned off to ensure we get the same results upon re-analysis (as Rsubread changes quite regularly, and I don't won't to have to keep on checking it). 
However, we still need to check that the code is correct. 
To do that, we run it with alignment turned back on.

```
cd aligntest
cat ../markdown/workflow.Rmd | sed "s/remap <- FALSE/remap <- TRUE/" > workflow.Rmd
echo "knitr::knit('workflow.Rmd')" | R --no-save --vanilla
cd -
```

The next step is to convert the markdown into a format of choice. 
Here, we convert it to HTML using render(), but it's also possible to convert to TeX via pandoc. 
The equivalent F1000 article was generated using the latter approach.

```
cd converted
cp ../markdown/workflow.md .
rsync -azv --delete ../markdown/figure . --exclude=".svn"
R -e "rmarkdown::render('workflow.md')"
cd -
```

We also need to convert the R markdown into something that BioC can compile.
This is done using the 'toBioC.sh' script, which flips some flags and modifies the YAML metadata so that it is compatible with BioC's build system.

```
cd bioC
./toBioC.sh
cd -
```
