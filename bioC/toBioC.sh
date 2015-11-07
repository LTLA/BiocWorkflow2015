# BioC's build system is more restrictive. We need to ensure that the
# 'biocbuild' flag is set to TRUE, to avoid the alignment steps that BioC can't
# do (as the command line utilities aren't available). Instead, they're
# starting from the BAM files. We also need to set keep.file to FALSE, in 
# order to tell the build system to re-download the files; and clear.memory 
# to TRUE, in order to get it to delete objects to save memory.

refdir='..'
cat ${refdir}/markdown/workflow.Rmd | sed "s/redownload <- FALSE/redownload <- TRUE/" | \
    sed "s/keep.files <- TRUE/keep.files <- FALSE/" | sed "s/clear.memory <- FALSE/clear.memory <- TRUE/" > temp.Rmd

# We also need to reformat the author listing and affiliations, because BioC's
# rendering system assumes that the author YAML field is a simple string.

allauthors=`cat temp.Rmd | grep "^  \- name:" | sed "s/.*: //g" | awk -vORS=", " '{ print $0 }' | sed 's/, $/\n/'`
allaffiliations=`cat temp.Rmd | grep "^    affiliation:" | sed "s/.*: //g" | sed "s/; /\n/g"  | awk '!a[$0]++' | awk -vORS="; " '{ print $0 }' | sed 's/; $/\n/'`

cat temp.Rmd | sed "s/^author:/author: $allauthors\nauthor_affiliation: $allaffiliations/" | grep -v "^  \- name:" | grep -v "^    affiliation" > temp2.Rmd

mv temp2.Rmd temp.Rmd

# We need to add a image size limiter to the HTML right after the YAML
# block, otherwise the image size on the page depends on its actual size.

cat temp.Rmd | perl -0777 -pe 's/---\n\n/---\n\n\n<style>\npre, img {\n  max-width: 100%;\n  display: block;\n}\n<\/style>\n\n\n\n/' > temp2.Rmd

mv temp2.Rmd temp.Rmd

# Deleting HTML comment tags, because they're not playing nice with
# BioC's markdown to HTML conversion.

cat temp.Rmd | perl -0777 -pe 's/<!--((?!-->).*)-->//s' > temp2.Rmd

mv temp2.Rmd temp.Rmd

# Moving to the BioC folder, if it isn't there already.

mv temp.Rmd ${refdir}/bioC/markdown/chipseq_db.Rmd
cp ${refdir}/converted/ref.bib ${refdir}/bioC/markdown 

exit 0

# This next section describes the BAM upload strategy to BioC's S3 servers.
# This is necessary in order for BioC's build system to access the files,
# as we can't unpack and align from the SRA files themselves.

cd ${refdir}/markdown
for x in `ls | egrep "\.(bam|bai)$"`
do
    echo $x
    aws s3 cp $x s3://chipseqdb-bamfiles --acl public-read 
done

aws s3 cp mm9-blacklist.bed s3://chipseqdb-bamfiles --acl public-read 
aws s3 cp mm9ToMm10.over.chain s3://chipseqdb-bamfiles --acl public-read 

