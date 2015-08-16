args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
subject <- args[1]

# Trying to run it.
rmd <- paste0(subject, ".Rmd")
final <- paste0(subject, ".md")
okay <- try(knitr::knit(rmd, final))

# Checking if it was successful.
if (!methods::is(okay, "try-error")) { q("no") }

# Otherwise, saving the current image.
cur.data <- paste0(subject, ".redone.Rda")
save.image(file=cur.data)

# Running a diff and identifying the final breakpoint (killing newlines added in *.md)
dummy.final <- paste0(final, ".tmp")
system(sprintf("sed -e :a -e '/^\\n*$/{$d;N;};/\\n$/ba' %s > %s", final, dummy.final))
x <- system(sprintf("diff %s %s | grep '^[1-9]' | tail -1", dummy.final, rmd), intern=TRUE)
keep <- as.integer(sub(",.*", "", sub(".*[cd]", "", x)))
unlink(dummy.final)

# Getting all lines.
all.lines <- readLines(rmd)
yet.run <- all.lines[keep:length(all.lines)]

all.requires <- all.lines[1:(keep-1)]
required <- grep("^require", all.requires)
all.requires <- all.requires[required]

# Constructing a new file with appropriate headers and 'require' calls.
if (grepl("_redo[0-9]+$", subject)) {
    curnum <- sub(".*_redo([0-9]+)$", "\\1", subject)
    curnum <- as.integer(curnum) + 1L
    newfile <- sub("_redo[0-9]+$", paste0("_redo", curnum), subject)
} else {
    newfile <- paste0(subject, "_redo1")
}

new.file <- paste0(newfile, ".Rmd")
writeLines(c(sprintf("
```{r, echo=FAlSE, results='hide', message=FALSE}
load('%s')
%s
opts_chunk$set(error=FALSE)
```
", cur.data, paste(all.requires, collapse="\n")),
    yet.run), con=new.file)

