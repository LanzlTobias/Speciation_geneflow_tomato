## Modified version of FD_escaneo_por_ventanas_union_de_ventanas_graficarlas.R
# available at: https://github.com/ericgonzalezs/Characterization_of_introgression_from_Zea_mays_ssp._mexicana_to_Mexican_highland_maize/blob/master/Introgression_analyses/
#
# For help run Rscript Fd_50_site_windows_v2.R --help

## Loading packages

if(require("optparse")){
  print("optparse is loaded correctly")
} else {
  print("trying to install optparse")
  install.packages("optparse")
  if(require("optparse")){
    print("optparse installed and loaded")
  } else {
    stop("could not install optparse")
  }
}

if(require("parallel")){
  print("parallel is loaded correctly")
} else {
  print("trying to install parallel")
  install.packages("parallel")
  if(require("parallel")){
    print("parallel installed and loaded")
  } else {
    stop("could not install parallel")
  }
}

if(require("dplyr")){
  print("dplyr is loaded correctly")
} else {
  print("trying to install dplyr")
  install.packages("dplyr")
  if(require("dplyr")){
    print("dplyr installed and loaded")
  } else {
    stop("could not install dplyr")
  }
}

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

## Define arguments

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="Input file with the ABBA BABA counts per sites created with ABBA_BABA.v1.pl during the ABBA_BABA-pipeline"),
  make_option(c("-o", "--output"), action="store", default="output", type='character',
              help="Output prefix, Default='output'"),
  make_option(c("-n", "--number_sites"), action="store", default=50, type='integer',
              help="Number of informative ABBA BABA sites per window, Default=50"),
  make_option(c("-q", "--quantile"), action="store", default=0.99, type='float',
              help="Quantile that should be taken for the unified windows (0.99 for the top 1 %), Default=0.99"),
  make_option(c("-t", "--threads"), action="store", default=1, type='integer',
              help="Number of threads that should be used, Default=1"))

## Reading in arguments

opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$i)) {
  stop("No input detected in --input")
}

## Loading input
data <- fread(opt$i,header=T,sep = "\t",fill=T)
data <- data[-c(nrow(data),(nrow(data)-1),(nrow(data)-2)),]

n_sites <- opt$n

data <- data[,c(1,2,15:22)]

Chr <- unique(data$CHR)

## Generating windows

generate_windows <- function(Chr,n_sites,data){
  df <- subset(data,data$CHR==Chr)
  table1 <- c()
  a <- 1
  for (i in 1:round((nrow(df)/n_sites))) {
    window <- df[a:c(a+(n_sites-1)),]
    a <- a+n_sites
    FdNum <- sum(window$FdNum)
    FdDenom <- sum(window$FdDenom)
    DNum <- sum(window$DNum)
    DDenom <- sum(window$DDenom)
    FhomNum <- sum(window$FhomNum)
    FhomDenom <- sum(window$FhomDenom)
    D <- DNum/DDenom
    FdRes <- FdNum/FdDenom
    window_start <- window[1,]
    window_end <- window[n_sites,]
    window_mid <- (window_start[,2]+window_end[,2])/2
    table1 <- rbind(table1,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
  }
  names(table1) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table1)[6:12],"Fd")
  table1
}

cl <- makeCluster(opt$t)
chr_list <- clusterApply(cl,
                         Chr,
                   generate_windows,
                   n_sites,
                   data)
stopCluster(cl)

table_without_filters <- rbindlist(chr_list)
table_with_filters <- rbindlist(chr_list) %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)

fwrite(table_without_filters, file = paste0(opt$o,"_50_site_WindowsWithoutFiltering.csv"))
fwrite(table_with_filters, file = paste0(opt$o,"_50_site_WindowsWithFiltering.csv"))

## Unification of windows

quant <- opt$q

unify_windows <- function(Chr,data,quant){
  table1 <- data %>%
    filter(Chromosome == Chr)
  threshold <- quantile(table1$Fd, quant)

  table1$YesOrNo <- ifelse(table1$Fd >= threshold, "Y","N")

  
  YesOrNo <- paste(table1$YesOrNo, collapse = "")
  starts <- gregexpr("NY",YesOrNo)
  stops <- gregexpr("YN",YesOrNo)
  
  window_start <- c()
  for (i in starts[[1]]){
    window_start <- append(window_start, table1$Window_start[i+1])
  }
  
  window_end <- c()
  for (i in stops[[1]]){
    window_end <- append(window_end, table1$Window_end[i])
  }
  if (length(window_start) > length(window_end)) {
    window_end[length(window_start)] <- tail(table1$Window_end, 1L)
  }
  if (length(window_start) < length(window_end)) {
    window_start <- c(head(table1$Window_start, 1L),window_start)
  }
  unified_windows <- cbind.data.frame(window_start,window_end)
  unified_windows$Chromosome <- paste0("Chr",chr_number)
  unified_windows$Window_size <- unified_windows$window_end - unified_windows$window_start
  unified_windows
}

cl <- makeCluster(opt$t)
uni_windows_list <- clusterApply(cl,
                                 Chr,
                           unify_windows,
                           table_with_filters,
                           quant)
stopCluster(cl)

uni_windows <- rbindlist(uni_windows_list)

fwrite(uni_windows, file = paste0(opt$o,"_50SNP_UnifiedWindows_",(1-quant),"percent.csv"))

