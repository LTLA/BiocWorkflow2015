# BioC's build system is more restrictive. We need to ensure that the
# 'biocbuild' flag is set to TRUE, to avoid the alignment steps that BioC can't
# do (as the command line utilities aren't available). Instead, they're
# starting from the BAM files. We also need to set keep.file to FALSE, in 
# order to tell the build system to re-download the files; and clear.memory 
# to TRUE, in order to get it to delete objects to save memory.

refdir='.'
cat ${refdir}/workflow.Rmd | sed "s/redownload <- FALSE/redownload <- TRUE/" | \
    sed "s/keep.files <- TRUE/keep.files <- FALSE/" | sed "s/clear.memory <- FALSE/clear.memory <- TRUE/" > temp.Rmd

# We need to add a image size limiter to the HTML right after the YAML
# block, otherwise the image size on the page depends on its actual size.

cat temp.Rmd | perl -0777 -pe 's/---\n\n/---\n\n\n<style>\npre, img {\n  max-width: 100%;\n  display: block;\n}\n<\/style>\n\n\n\n/' > temp2.Rmd

mv temp2.Rmd temp.Rmd

# Deleting HTML comment tags, because they're not playing nice with
# BioC's markdown to HTML conversion.

cat temp.Rmd | perl -0777 -pe 's/<!--((?!-->).*)-->//s' > temp2.Rmd

mv temp2.Rmd temp.Rmd

# Moving to the BioC folder, if it isn't there already.

reldir='../bioc'
mv temp.Rmd ${reldir}/chipseq_db.Rmd
cp ${refdir}/ref.bib ${reldir}

exit 0

# This next section describes the BAM upload strategy to BioC's S3 servers.
# This is necessary in order for BioC's build system to access the files,
# as we can't unpack and align from the SRA files themselves.

for x in `ls | egrep "\.(bam|bai)$"`
do
    echo $x
    aws s3 cp $x s3://chipseqdb-bamfiles --acl public-read 
done

aws s3 cp mm10.blacklist.bed.gz s3://chipseqdb-bamfiles --acl public-read 

