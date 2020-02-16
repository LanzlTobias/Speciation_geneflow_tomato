#########################
#Input
#########################
#First the unzipped vcf
#Second is the gff3 file
#Third is the fasta file of the reference genome
#Fourth is a file containing all the population names (one population per line)

## Loading required packages
if(require("data.table")){
  print("data.table is loaded correctly")
} else {
  print("trying to install data.table")
  install.packages("data.table")
  if(require("data.table")){
    print("data.table installed and loaded")
  } else {
    stop("could not install data.table")
  }
}
if(require("tidyverse")){
  print("tidyverse is loaded correctly")
} else {
  print("trying to install tidyverse")
  install.packages("tidyverse")
  if(require("tidyverse")){
    print("tidyverse installed and loaded")
  } else {
    stop("could not install tidyverse")
  }
}

if(require("BiocManager")){
  print("BiocManager is loaded correctly")
} else {
  print("trying to install BiocManager")
  install.packages("BiocManager")
  if(require("BiocManager")){
    print("BiocManager installed and loaded")
  } else {
    stop("could not install BiocManager")
  }
}

if(require("VariantAnnotation")){
  print("VariantAnnotation is loaded correctly")
} else {
  print("trying to install VariantAnnotation")
  BiocManager::install("VariantAnnotation")
  if(require("VariantAnnotation")){
    print("VariantAnnotation installed and loaded")
  } else {
    stop("could not install VariantAnnotation")
  }
}

if(require("GenomicFeatures")){
  print("GenomicFeatures is loaded correctly")
} else {
  print("trying to install GenomicFeatures")
  BiocManager::install("GenomicFeatures")
  if(require("GenomicFeatures")){
    print("GenomicFeatures installed and loaded")
  } else {
    stop("could not install GenomicFeatures")
  }
}

if(require("Rsamtools")){
  print("Rsamtools is loaded correctly")
} else {
  print("trying to install Rsamtools")
  BiocManager::install("Rsamtools")
  if(require("Rsamtools")){
    print("Rsamtools installed and loaded")
  } else {
    stop("could not install Rsamtools")
  }
}



## Reading in and preparing data
input <- commandArgs(trailingOnly = TRUE)


fa <- open(FaFile(input[3]))
txdb <- makeTxDbFromGFF(input[2],format="gff3")
vcf_object <- readVcf(file=input[1])
effects <- predictCoding(vcf_object,txdb,fa)

id_from_effects <- names(ranges(effects))
id_from_vcf <- rownames(info(vcf_object))
m <- match(id_from_vcf,id_from_effects)

info(vcf_object)$CSQ <- effects$CONSEQUENCE[m]
vcf_objects <- vcf_object[!is.na(m)]

x <- effects$GENEID[m]
x[is.na(x)] <- "none"

names(vcf_object) <- x

read.vcf = function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}
my_data <- read.vcf(writeVcf(vcf_object, paste0("wAnnot_",input[1])),header=TRUE,stringsAsFactors=FALSE)

my_data$CSQ <- effects$CONSEQUENCE[m]


### Testing all loci
population <- scan(input[4],character())


mk_test <- function(pop){
df <- my_data[,c(1:3,ncol(my_data),grep(pop,colnames(my_data)))] %>%
  filter(ID != "none")

df <- df %>%
  pivot_longer(cols=5:ncol(df),names_to = "ind",values_to = "genotype")
df$genotype <- gsub(":.+","",df$genotype)
df$genotype <- gsub("\\|","/",df$genotype)
df <- df %>%
  separate(genotype,into=c("A1","A2"),sep="/") %>%
  filter(!is.na(A2)) %>%
  pivot_longer(cols=A1:A2, names_to = "Allele",values_to = "count")

result <- df %>%
  group_by(CHROM,POS,ID,CSQ) %>%
  summarize(N=sum(as.numeric(count))) %>%
  ungroup() %>%
  filter(N > 0) %>%
  mutate(state=ifelse(N == 10, "fixed","polymorphic")) %>%
  group_by(ID) %>%
  summarize(Fn = sum(state == "fixed" & CSQ == "nonsynonymous"),
            Fs = sum(state == "fixed" & CSQ == "synonymous"),
            Pn = sum(state == "polymorphic" & CSQ == "nonsynonymous"),
            Ps = sum(state == "polymorphic" & CSQ == "synonymous"),
            between=Fn/Fs, within = Pn/Ps,
            NI=within/between,
            nsyn = Fn/Pn, syn=Fs/Ps,
            neutral = nsyn/syn,
            DoS = Fn/(Fn+Fs)-Pn/(Pn+Ps))

fisher_test <- function(n) {
  fn <- result$Fn[n]
  fs <- result$Fs[n]
  pn <- result$Pn[n]
  ps <- result$Ps[n]

  x <- matrix(c(fn,pn,fs,ps),nrow=2,ncol=2,byrow=TRUE)
  fisher.test(x)$p.value
}
result$p_value <- sapply(1:nrow(result),fisher_test)

fwrite(result,file=paste0(pop,"_",str_remove(input[1],".vcf"),"_allLoci.csv"))
result
}
mk_results <- lapply(population,mk_test)
names(mk_results) <- population
mk_df <- rbindlist(mk_results,idcol="pop")
fwrite(mk_df,file=paste0("MK_",str_remove(input[1],".vcf"),"_allLoci.csv"))
