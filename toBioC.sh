# BioC's build system is more restrictive. We need to flip the 'biocbuild' flag 
# to TRUE, to avoid the alignment steps that BioC can't do (as the command line 
# utilities aren't available). Instead, they're starting from the BAM files.

cat markdown/workflow.Rmd | sed "s/biocbuild <- FALSE/biocbuild <- TRUE/" > temp.Rmd

# We also need to reformat the author listing and affiliations, because BioC's
# rendering system assumes that the author YAML field is a simple string.

allauthors=`cat temp.Rmd | grep "^  \- name:" | sed "s/.*: //g" | awk -vORS=", " '{ print $0 }' | sed 's/, $/\n/'`
allaffiliations=`cat temp.Rmd | grep "^    affiliation:" | sed "s/.*: //g" | sed "s/; /\n/g"  | awk '!a[$0]++' | awk -vORS="; " '{ print $0 }' | sed 's/; $/\n/'`

cat temp.Rmd | sed "s/^author:/author: $allauthors\nauthor_affiliation: $allaffiliations/" | grep -v "^  \- name:" | grep -v "^    affiliation" > temp2.Rmd

mv temp2.Rmd temp.Rmd

# We need to add a image size limiter to the HTML right after the YAML
# block, otherwise the image size on the page depends on its actual size.

cat temp.Rmd | sed -e '2,$s/---/---\n\n\n<style>\npre, img {\n  max-width: 100%;\n  display: block;\n}\n<\/style>\n\n/' > temp2.Rmd

mv temp2.Rmd temp.Rmd

# Moving to the BioC folder, if it isn't there already.

mv temp.Rmd bioC/chipseq_db.Rmd
cp converted/ref.bib bioC

exit 0

# This next section describes the BAM upload strategy to BioC's S3 servers.
# This is necessary in order for BioC's build system to access the files,
# as we can't unpack and align from the SRA files themselves.

for x in `ls markdown/* | egrep "\.(bam|bai)$"`
do
    echo $x
    aws s3 cp $x s3://chipseqdb-bamfiles --acl public-read 
done

aws s3 cp markdown/mm9-blacklist.bed s3://chipseqdb-bamfiles --acl public-read 
aws s3 cp markdown/mm9ToMm10.over.chain s3://chipseqdb-bamfiles --acl public-read 

