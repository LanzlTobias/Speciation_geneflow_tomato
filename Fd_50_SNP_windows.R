## Modified version of FD_escaneo_por_ventanas_union_de_ventanas_graficarlas.R
# available at: https://github.com/ericgonzalezs/Characterization_of_introgression_from_Zea_mays_ssp._mexicana_to_Mexican_highland_maize/blob/master/Introgression_analyses/
#
#Arguments for the script:
#First argument is the file with the ABBA BABA counts per sites created with ABBA_BABA.v1.pl during the ABBA_BABA-pipeline
#Second argument is the output prefix

#Loading packages
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


## Loading input
inputfile <- commandArgs(trailingOnly = TRUE)
data <- fread(inputfile[1],header=T,sep = "\t",fill=T)
data <- data[-c(nrow(data),(nrow(data)-1),(nrow(data)-2)),]


Fs <- data[,c(1,2,15:22)]
Fs_1 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[1])
Fs_2 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[2])
Fs_3 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[3])
Fs_4 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[4])
Fs_5 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[5])
Fs_6 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[6])
Fs_7 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[7])
Fs_8 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[8])
Fs_9 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[9])
Fs_10 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[10])
Fs_11 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[11])
Fs_12 <- subset(Fs, Fs$CHR==unique(Fs$CHR)[12])


#####
#50 SNP windows

#chr1
table1 <- c()
a <- 1
for (i in 1:round((nrow(Fs_1)/50))) {
  window <- Fs_1[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table1 <- rbind(table1,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table1) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table1)[6:12],"Fd")
#chr2
table2 <- c()
a <- 1
for (i in 1:round((nrow(Fs_2)/50))) {
  window <- Fs_2[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table2 <- rbind(table2,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table2) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table2)[6:12],"Fd")

#chr3
table3 <- c()
a <- 1
for (i in 1:round((nrow(Fs_3)/50))) {
  window <- Fs_3[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table3 <- rbind(table3,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table3) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table3)[6:12],"Fd")

#chr4
table4 <- c()
a <- 1
for (i in 1:round((nrow(Fs_4)/50))) {
  window <- Fs_4[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table4 <- rbind(table4,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table4) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table4)[6:12],"Fd")

#chr5
table5 <- c()
a <- 1
for (i in 1:round((nrow(Fs_5)/50))) {
  window <- Fs_5[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table5 <- rbind(table5,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table5) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table5)[6:12],"Fd")


#chr6
table6 <- c()
a <- 1
for (i in 1:round((nrow(Fs_6)/50))) {
  window <- Fs_6[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table6 <- rbind(table6,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table6) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table6)[6:12],"Fd")

#chr7
table7 <- c()
a <- 1
for (i in 1:round((nrow(Fs_7)/50))) {
  window <- Fs_7[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table7 <- rbind(table7,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table7) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table7)[6:12],"Fd")

#chr8
table8 <- c()
a <- 1
for (i in 1:round((nrow(Fs_8)/50))) {
  window <- Fs_8[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table8 <- rbind(table8,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table8) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table8)[6:12],"Fd")

#chr9
table9 <- c()
a <- 1
for (i in 1:round((nrow(Fs_9)/50))) {
  window <- Fs_9[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table9 <- rbind(table9,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table9) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table9)[6:12],"Fd")

#chr10
table10 <- c()
a <- 1
for (i in 1:round((nrow(Fs_10)/50))) {
  window <- Fs_10[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table10 <- rbind(table10,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table10) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table10)[6:12],"Fd")

#chr11
table11 <- c()
a <- 1
for (i in 1:round((nrow(Fs_11)/50))) {
  window <- Fs_11[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table11 <- rbind(table11,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table11) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table11)[6:12],"Fd")

#chr12
table12 <- c()
a <- 1
for (i in 1:round((nrow(Fs_12)/50))) {
  window <- Fs_12[a:c(a+49),]
  a <- a+50
  FdNum <- sum(window$FdNum)
  FdDenom <- sum(window$FdDenom)
  DNum <- sum(window$DNum)
  DDenom <- sum(window$DDenom)
  FhomNum <- sum(window$FhomNum)
  FhomDenom <- sum(window$FhomDenom)
  D <- DNum/DDenom
  FdRes <- FdNum/FdDenom
  window_start <- window[1,]
  window_end <- window[50,]
  window_mid <- (window_start[,2]+window_end[,2])/2
  table12 <- rbind(table12,data.frame(window[1,1],window_start[,2],window_end[,2],window_end[,2]-window_start[,2],window_mid,FdNum, FdDenom, DNum, DDenom, FhomNum, FhomDenom, D, FdRes))
}
names(table12) <- c("Chromosome","Window_start","Window_end","Window_size","Window_mid",names(table12)[6:12],"Fd")

#####
#Filter for D < 0 and Fd between 0 and 1
table1_subset <- table1 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table2_subset <- table2 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table3_subset <- table3 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table4_subset <- table4 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table5_subset <- table5 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table6_subset <- table6 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table7_subset <- table7 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table8_subset <- table8 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table9_subset <- table9 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table10_subset <- table10 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table11_subset <- table11 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)
table12_subset <- table12 %>%
  filter(D > 0, Fd >= 0 & Fd <= 1)

table_without_filters <- rbind.data.frame(table1, table2, table3, table4, table5, table6, table7, table8, table9, table10,table11,table12)
table_with_filters <- rbind.data.frame(table1_subset, table2_subset, table3_subset, table4_subset, table5_subset, table6_subset, table7_subset, table8_subset, table9_subset, table10_subset,table11_subset,table12_subset)

fwrite(table_without_filters, file = paste0(inputfile[2],"_50SNP_WindowsWithoutFiltering.csv"))
fwrite(table_with_filters, file = paste0(inputfile[2],"_50SNP_WindowsWithFiltering.csv"))

###################################################
#write.csv(table_with_filters, file = "")
##################################################

#####
#Unification of windows

#Chr1
quant_90 <- quantile(table1_subset$Fd, 0.90)
quant_99 <- quantile(table1_subset$Fd, 0.99)

table1_subset$YesOrNo_90 <- ifelse(table1_subset$Fd >= quant_90, "Y","N")
table1_subset$YesOrNo_99 <- ifelse(table1_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table1_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table1_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table1_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table1_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table1_subset$Window_start, 1L),window_start)
}
Chr1_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr1_unified_windows_90$Chromosome <- "Chr1"
Chr1_unified_windows_90$Window_size <- Chr1_unified_windows_90$window_end - Chr1_unified_windows_90$window_start


YesOrNo <- paste(table1_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table1_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table1_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table1_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table1_subset$Window_start, 1L),window_start)
}
Chr1_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr1_unified_windows_99$Chromosome <- "Chr1"
Chr1_unified_windows_99$Window_size <- Chr1_unified_windows_99$window_end - Chr1_unified_windows_99$window_start

#Chr2
quant_90 <- quantile(table2_subset$Fd, 0.90)
quant_99 <- quantile(table2_subset$Fd, 0.99)

table2_subset$YesOrNo_90 <- ifelse(table2_subset$Fd >= quant_90, "Y","N")
table2_subset$YesOrNo_99 <- ifelse(table2_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table2_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table2_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table2_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table2_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table2_subset$Window_start, 1L),window_start)
}
Chr2_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr2_unified_windows_90$Chromosome <- "Chr2"
Chr2_unified_windows_90$Window_size <- Chr2_unified_windows_90$window_end - Chr2_unified_windows_90$window_start


YesOrNo <- paste(table2_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table2_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table2_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table2_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table2_subset$Window_start, 1L),window_start)
}
Chr2_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr2_unified_windows_99$Chromosome <- "Chr2"
Chr2_unified_windows_99$Window_size <- Chr2_unified_windows_99$window_end - Chr2_unified_windows_99$window_start

#Chr3
quant_90 <- quantile(table3_subset$Fd, 0.90)
quant_99 <- quantile(table3_subset$Fd, 0.99)

table3_subset$YesOrNo_90 <- ifelse(table3_subset$Fd >= quant_90, "Y","N")
table3_subset$YesOrNo_99 <- ifelse(table3_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table3_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table3_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table3_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table3_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table3_subset$Window_start, 1L),window_start)
}
Chr3_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr3_unified_windows_90$Chromosome <- "Chr3"
Chr3_unified_windows_90$Window_size <- Chr3_unified_windows_90$window_end - Chr3_unified_windows_90$window_start


YesOrNo <- paste(table3_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table3_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table3_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table3_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table3_subset$Window_start, 1L),window_start)
}
Chr3_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr3_unified_windows_99$Chromosome <- "Chr3"
Chr3_unified_windows_99$Window_size <- Chr3_unified_windows_99$window_end - Chr3_unified_windows_99$window_start

#Chr4
quant_90 <- quantile(table4_subset$Fd, 0.90)
quant_99 <- quantile(table4_subset$Fd, 0.99)

table4_subset$YesOrNo_90 <- ifelse(table4_subset$Fd >= quant_90, "Y","N")
table4_subset$YesOrNo_99 <- ifelse(table4_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table4_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table4_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table4_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table4_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table4_subset$Window_start, 1L),window_start)
}
Chr4_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr4_unified_windows_90$Chromosome <- "Chr4"
Chr4_unified_windows_90$Window_size <- Chr4_unified_windows_90$window_end - Chr4_unified_windows_90$window_start


YesOrNo <- paste(table4_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table4_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table4_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table4_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table4_subset$Window_start, 1L),window_start)
}
Chr4_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr4_unified_windows_99$Chromosome <- "Chr4"
Chr4_unified_windows_99$Window_size <- Chr4_unified_windows_99$window_end - Chr4_unified_windows_99$window_start

#Chr5
quant_90 <- quantile(table5_subset$Fd, 0.90)
quant_99 <- quantile(table5_subset$Fd, 0.99)

table5_subset$YesOrNo_90 <- ifelse(table5_subset$Fd >= quant_90, "Y","N")
table5_subset$YesOrNo_99 <- ifelse(table5_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table5_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table5_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table5_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table5_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table5_subset$Window_start, 1L),window_start)
}
Chr5_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr5_unified_windows_90$Chromosome <- "Chr5"
Chr5_unified_windows_90$Window_size <- Chr5_unified_windows_90$window_end - Chr5_unified_windows_90$window_start


YesOrNo <- paste(table5_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table5_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table5_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table5_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table5_subset$Window_start, 1L),window_start)
}
Chr5_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr5_unified_windows_99$Chromosome <- "Chr5"
Chr5_unified_windows_99$Window_size <- Chr5_unified_windows_99$window_end - Chr5_unified_windows_99$window_start

#Chr6
quant_90 <- quantile(table6_subset$Fd, 0.90)
quant_99 <- quantile(table6_subset$Fd, 0.99)

table6_subset$YesOrNo_90 <- ifelse(table6_subset$Fd >= quant_90, "Y","N")
table6_subset$YesOrNo_99 <- ifelse(table6_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table6_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table6_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table6_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table6_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table6_subset$Window_start, 1L),window_start)
}
Chr6_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr6_unified_windows_90$Chromosome <- "Chr6"
Chr6_unified_windows_90$Window_size <- Chr6_unified_windows_90$window_end - Chr6_unified_windows_90$window_start


YesOrNo <- paste(table6_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table6_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table6_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table6_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table6_subset$Window_start, 1L),window_start)
}
Chr6_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr6_unified_windows_99$Chromosome <- "Chr6"
Chr6_unified_windows_99$Window_size <- Chr6_unified_windows_99$window_end - Chr6_unified_windows_99$window_start

#Chr7
quant_90 <- quantile(table7_subset$Fd, 0.90)
quant_99 <- quantile(table7_subset$Fd, 0.99)

table7_subset$YesOrNo_90 <- ifelse(table7_subset$Fd >= quant_90, "Y","N")
table7_subset$YesOrNo_99 <- ifelse(table7_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table7_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table7_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table7_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table7_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table7_subset$Window_start, 1L),window_start)
}
Chr7_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr7_unified_windows_90$Chromosome <- "Chr7"
Chr7_unified_windows_90$Window_size <- Chr7_unified_windows_90$window_end - Chr7_unified_windows_90$window_start


YesOrNo <- paste(table7_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table7_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table7_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table7_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table7_subset$Window_start, 1L),window_start)
}
Chr7_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr7_unified_windows_99$Chromosome <- "Chr7"
Chr7_unified_windows_99$Window_size <- Chr7_unified_windows_99$window_end - Chr7_unified_windows_99$window_start

#Chr8
quant_90 <- quantile(table8_subset$Fd, 0.90)
quant_99 <- quantile(table8_subset$Fd, 0.99)

table8_subset$YesOrNo_90 <- ifelse(table8_subset$Fd >= quant_90, "Y","N")
table8_subset$YesOrNo_99 <- ifelse(table8_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table8_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table8_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table8_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table8_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table8_subset$Window_start, 1L),window_start)
}
Chr8_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr8_unified_windows_90$Chromosome <- "Chr8"
Chr8_unified_windows_90$Window_size <- Chr8_unified_windows_90$window_end - Chr8_unified_windows_90$window_start


YesOrNo <- paste(table8_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table8_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table8_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table8_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table8_subset$Window_start, 1L),window_start)
}
Chr8_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr8_unified_windows_99$Chromosome <- "Chr8"
Chr8_unified_windows_99$Window_size <- Chr8_unified_windows_99$window_end - Chr8_unified_windows_99$window_start

#Chr9
quant_90 <- quantile(table9_subset$Fd, 0.90)
quant_99 <- quantile(table9_subset$Fd, 0.99)

table9_subset$YesOrNo_90 <- ifelse(table9_subset$Fd >= quant_90, "Y","N")
table9_subset$YesOrNo_99 <- ifelse(table9_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table9_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)


window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table9_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table9_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table9_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table9_subset$Window_start, 1L),window_start)
}
Chr9_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr9_unified_windows_90$Chromosome <- "Chr9"
Chr9_unified_windows_90$Window_size <- Chr9_unified_windows_90$window_end - Chr9_unified_windows_90$window_start


YesOrNo <- paste(table9_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table9_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table9_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table9_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table9_subset$Window_start, 1L),window_start)
}
Chr9_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr9_unified_windows_99$Chromosome <- "Chr9"
Chr9_unified_windows_99$Window_size <- Chr9_unified_windows_99$window_end - Chr9_unified_windows_99$window_start

#Chr10
quant_90 <- quantile(table10_subset$Fd, 0.90)
quant_99 <- quantile(table10_subset$Fd, 0.99)

table10_subset$YesOrNo_90 <- ifelse(table10_subset$Fd >= quant_90, "Y","N")
table10_subset$YesOrNo_99 <- ifelse(table10_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table10_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table10_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table10_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table10_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table10_subset$Window_start, 1L),window_start)
}
Chr10_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr10_unified_windows_90$Chromosome <- "Chr10"
Chr10_unified_windows_90$Window_size <- Chr10_unified_windows_90$window_end - Chr10_unified_windows_90$window_start


YesOrNo <- paste(table10_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table10_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table10_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table10_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table10_subset$Window_start, 1L),window_start)
}
Chr10_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr10_unified_windows_99$Chromosome <- "Chr10"
Chr10_unified_windows_99$Window_size <- Chr10_unified_windows_99$window_end - Chr10_unified_windows_99$window_start

#Chr11
quant_90 <- quantile(table11_subset$Fd, 0.90)
quant_99 <- quantile(table11_subset$Fd, 0.99)

table11_subset$YesOrNo_90 <- ifelse(table11_subset$Fd >= quant_90, "Y","N")
table11_subset$YesOrNo_99 <- ifelse(table11_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table11_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table11_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table11_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table11_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table11_subset$Window_start, 1L),window_start)
}
Chr11_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr11_unified_windows_90$Chromosome <- "Chr11"
Chr11_unified_windows_90$Window_size <- Chr11_unified_windows_90$window_end - Chr11_unified_windows_90$window_start


YesOrNo <- paste(table11_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table11_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table11_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table11_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table11_subset$Window_start, 1L),window_start)
}
Chr11_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr11_unified_windows_99$Chromosome <- "Chr11"
Chr11_unified_windows_99$Window_size <- Chr11_unified_windows_99$window_end - Chr11_unified_windows_99$window_start

#Chr12
quant_90 <- quantile(table12_subset$Fd, 0.90)
quant_99 <- quantile(table12_subset$Fd, 0.99)

table12_subset$YesOrNo_90 <- ifelse(table12_subset$Fd >= quant_90, "Y","N")
table12_subset$YesOrNo_99 <- ifelse(table12_subset$Fd >= quant_99, "Y","N")


YesOrNo <- paste(table12_subset$YesOrNo_90, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table12_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table12_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table12_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table12_subset$Window_start, 1L),window_start)
}
Chr12_unified_windows_90 <- cbind.data.frame(window_start,window_end)
Chr12_unified_windows_90$Chromosome <- "Chr12"
Chr12_unified_windows_90$Window_size <- Chr12_unified_windows_90$window_end - Chr12_unified_windows_90$window_start


YesOrNo <- paste(table12_subset$YesOrNo_99, collapse = "")
starts <- gregexpr("NY",YesOrNo)
stops <- gregexpr("YN",YesOrNo)

window_start <- c()
for (i in starts[[1]]){
  window_start <- append(window_start, table12_subset$Window_start[i+1])
}

window_end <- c()
for (i in stops[[1]]){
  window_end <- append(window_end, table12_subset$Window_end[i])
}
if (length(window_start) > length(window_end)) {
  window_end[length(window_start)] <- tail(table12_subset$Window_end, 1L)
}
if (length(window_start) < length(window_end)) {
  window_start <- c(head(table12_subset$Window_start, 1L),window_start)
}
Chr12_unified_windows_99 <- cbind.data.frame(window_start,window_end)
Chr12_unified_windows_99$Chromosome <- "Chr12"
Chr12_unified_windows_99$Window_size <- Chr12_unified_windows_99$window_end - Chr12_unified_windows_99$window_start


Unified_windows_90 <- bind_rows(Chr1_unified_windows_90,Chr2_unified_windows_90,Chr3_unified_windows_90,Chr4_unified_windows_90,Chr5_unified_windows_90,Chr6_unified_windows_90,Chr7_unified_windows_90,Chr8_unified_windows_90,Chr9_unified_windows_90,Chr10_unified_windows_90,Chr11_unified_windows_90,Chr12_unified_windows_90)
Unified_windows_99 <- bind_rows(Chr1_unified_windows_99,Chr2_unified_windows_99,Chr3_unified_windows_99,Chr4_unified_windows_99,Chr5_unified_windows_99,Chr6_unified_windows_99,Chr7_unified_windows_99,Chr8_unified_windows_99,Chr9_unified_windows_99,Chr10_unified_windows_99,Chr11_unified_windows_99,Chr12_unified_windows_99)
fwrite(Unified_windows_90, file = paste0(inputfile[2],"_50SNP_UnifiedWindows_10percent.csv"))
fwrite(Unified_windows_99, file = paste0(inputfile[2],"_50SNP_UnifiedWindows_1percent.csv"))
