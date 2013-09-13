#!/bin/bash

snpfiles=$1


base=~/dwork/snpspec

Hu2011=(
$base/data/Hu2011/RA_SNPs.txt
$base/data/Hu2011/Crohns_SNPs.txt
$base/data/Hu2011/SLE_SNPs.txt
)

MarkLee=(
$base/data/MarkLee/AnkylosingSpondylitis.Lin2011.txt
$base/data/MarkLee/CeliacDisease.Dubois2010.txt
$base/data/MarkLee/Psoriasis.Tsoi2012.txt
$base/data/MarkLee/Leprosy.Zhang2009.txt
$base/data/MarkLee/MultipleSclerosis.Sawcer2011.txt
$base/data/MarkLee/Type1Diabetes.Barrett2009.txt
)

Niko=(
$base/data/NikoPatsopoulos/ai.0.1_ms-sp.0.1_effects.txt
$base/data/NikoPatsopoulos/ai.0.1.effects.txt
$base/data/NikoPatsopoulos/ms-sp.0.1.effects.txt
)

Towfique=(
$base/data/Towfique/AD_ALL.txt
$base/data/Towfique/alzheimers.txt
$base/data/Towfique/ADGC_EUR_GWAS.txt
$base/data/Towfique/AD_1KG_1e-6.txt
$base/data/Towfique/AD_1KG_1e-7.txt
)

JessicavanSetten=(
$base/data/JessicavanSetten/QTIGC_SNPs.txt
$base/data/JessicavanSetten/61_SNPs.txt
)

Random=(
$base/data/Random/*.txt
)

#data/BrainAtlas/schizophrenia_snps.snpedia.txt
#data/BrainAtlas/neuroticism_snps.Calboli2010.txt
snpfiles=(
    $base/data/Harst2012/RBC_SNPs.txt
    $base/data/IGAS2013/AS_SNPs.txt
    ${Hu2011[*]}
    ${MarkLee[*]}
    ${Niko[*]}
    ${Towfique[*]}
    ${JessicavanSetten[*]}
)

immgen2010=$base/data/ImmGen/ImmGen2010.gct.gz
immgen2012=$base/data/ImmGen/ImmGen2012.gct.gz
novartis=$base/data/Novartis/Novartis2011.gct.gz
go2006=$base/data/GO/GO2006/allHumanGO.gct.gz
go2006_500=$base/data/GO/GO2006/allHumanGO.gte500_and_lte1000.gct.gz
go2013=$base/data/GO/allHumanGO.gct.gz

#snps=$base/data/Hu2011/RA_SNPs.txt
#condition=$base/data/GO/first_1000.txt

allsnps=$(echo ${snpfiles[*]} | tr ' ' ',')

#for snps in ${snpfiles[*]} # ${Random[*]}
for snps in $allsnps
do
    #echo $snps
    for expression in $immgen2010 $immgen2012 $novartis $go2013
    do
        echo $(date) $expression

        #a=$(echo $snps | perl -ne 'm`([^/]+)/([^/]+)\.txt` && print "$1_$2"')
        out=out/batch/$(date +%Y%m%d) # /${a}_$(basename $expression .gct.gz)
        mkdir -p $out

        #    --condition $condition
        options=(
            --snps $snps
            --expression $expression
            --gene-intervals $base/reference/genes.hg19.ucsc.bed
            --null-snps $base/reference/TGP.pruned.txt.gz
            --out $out
            --snp-intervals $base/reference/TGP.bed.gz
            --slop 250e3
            --processes 6
            --min-observations 50
            --max-iterations 1e6
        )
        time ./snpspec ${options[*]} | tee $out/log.txt
    done
done
