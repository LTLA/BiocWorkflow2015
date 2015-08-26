pandoc -t latex workflow.md --bibliography ref.bib --biblatex \
    --template ~/Software/R/current/library/rmarkdown/rmd/latex/default.tex \
    --highlight-style tango > temp.tex

cat temp.tex | sed  "s/caption..textbf{Figure [^:]*:} /caption{/" > temp2.tex
mv temp2.tex temp.tex

cat temp.tex | sed "s/figure\///" > temp2.tex
mv temp2.tex temp.tex	

cat temp.tex | sed "s/{biblatex}/[style=authoryear]{biblatex}/" > temp2.tex
mv temp2.tex temp.tex

cat temp.tex | sed -n -e '/\\section{Introduction}/,$p' > temp2.tex
mv temp2.tex temp.tex

cat preamble.tex temp.tex > workflow.tex
rm temp.tex
