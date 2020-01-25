## Loading necessary packages
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

if(require("data.table")){
  print("data.table is loaded correctly")
} else {
  print("trying to install data.table")
  install.packages("data.table")
  if(require(data.table)){
    print("data.table installed and loaded")
  } else {
    stop("could not install data.table")
  }
}



popstat_ANGSD <- fread("popgenWindows.csv")
chr_old <- unique(popstat_ANGSD$scaffold)
chr_new <- c(paste0(rep("Chr",12),1:12))

for (i in 1:12){
  popstat_ANGSD$scaffold[popstat_ANGSD$scaffold == chr_old[i]] <- chr_new[i]
}
popstat_ANGSD$scaffold <- factor(popstat_ANGSD$scaffold, levels = chr_new)

positions <- grep("LA1963",colnames(popstat_ANGSD))[-1]
div_df <- popstat_ANGSD[,c(1:13,..positions)] %>%
  mutate(Dxy = dxy_LA1963_LA_Peruvianum,Fst = Fst_LA1963_LA_Peruvianum)
colnames(div_df)[1] <- "Chromosome"

Zdiv_df <- div_df %>%
  filter(!is.na(Fst)) %>%
  mutate(ZFst = (Fst-mean(Fst))/sd(Fst),
         ZDxy = (Dxy - mean(Dxy))/sd(Dxy),
         Z_pi_LA1963 = (pi_LA1963-mean(pi_LA1963))/sd(pi_LA1963),
         Z_pi_LA_Peruvianum = (pi_LA_Peruvianum-mean(pi_LA_Peruvianum,na.rm=TRUE))/sd(pi_LA_Peruvianum,na.rm=TRUE),
         position=cumsum(as.numeric(mid))) 

speciation_df_single <- Zdiv_df %>%
  filter(sites >= 100,
         ZDxy > 2,
         ZFst > 2 ) %>%
  mutate(region_start = ((start+end)/2)-250000,
         region_end = ((start+end)/2)+250000,
         Region_number = row_number())

filter_regions <- function(n,df,regions){
  x <- regions[n,]
  df %>%
    filter(Chromosome==x$Chromosome & mid >= x$start & mid <= x$end)
}


### SNP density filter

SNP_density<- fread("solanum_10kb.snpden")
chr_old <- unique(SNP_density$CHROM)
chr_new <- c(paste0("Chr",1:12))

for (i in 1:12){
  SNP_density$CHROM[SNP_density$CHROM == chr_old[i]] <- chr_new[i]
}
SNP_density$CHROM <- factor(SNP_density$CHROM, levels = chr_new)

SNP_density <- SNP_density %>%
  mutate(Chromosome = CHROM, mid = BIN_START + 5000)


snpden_s_list <- lapply(1:nrow(speciation_df_single),
                        filter_regions,
                        SNP_density,
                        speciation_df_single)

snpden_df <-     rbindlist(snpden_s_list,idcol="Region_number")

filter_df <- snpden_df %>%
  group_by(Region_number) %>%
  summarize(density = mean(`VARIANTS/KB`),n_10kb=n(),n_below=sum(`VARIANTS/KB` < 10)) %>%
  ungroup() %>%
  mutate(filter_step1 = n_below/n_10kb,filter_step2= filter_step1 >= 0.5) %>%
  filter(filter_step2)

### Preparation of 10 kb popstat file

popstat_10kb <- fread("popgenWindows_10Kb.csv")
chr_old <- unique(popstat_10kb$scaffold)
chr_new <- c(paste0(rep("Chr",12),1:12))
for (i in 1:12){
  popstat_10kb$scaffold[popstat_10kb$scaffold == chr_old[i]] <- chr_new[i]
}
popstat_10kb$scaffold <- factor(popstat_10kb$scaffold, levels = chr_new)

positions <- grep("LA1963",colnames(popstat_10kb))[-1]
div_10kb <- popstat_10kb[,c(1:13,..positions)]
colnames(div_10kb)[1] <- "Chromosome"

scale2 <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm=TRUE)
div_10kb_Z <- div_10kb %>%
  mutate_at(names(div_10kb)[14:27],scale2)

div_10kb_Z <- div_10kb_Z %>%
  gather(key="Variable",value="Value",-(Chromosome:sites)) %>%
  separate(col="Variable",into = c("Variable","Pi_pop","Population",NA)) %>%
  mutate(Out = Variable != "pi" & Pi_pop != "LA1963") %>%
  filter(!Out) %>%
  mutate(Population = ifelse(is.na(Population) | Population == "Peruvianum",Pi_pop,Population),
         Value = ifelse(Value < 0,0,Value)) %>%
  filter(sites >= 10)

### Speciation pattern filter

filter_regions <- function(n,df,regions){
  x <- regions[n,]
  df %>%
    filter(Chromosome==x$Chromosome & mid >= x$start & mid <= x$end)
}

s_regions_list <- lapply(1:nrow(speciation_df_single),
                         filter_regions,
                         div_10kb_Z,
                         speciation_df_single)


s_regions_df <- rbindlist(s_regions_list,idcol="Region_number")

pattern_summary <- filter(s_regions_df,!(Population %in% c("LA","PI")),                                                        !(Variable %in% c("Tajima_D","LD"))) %>%
  group_by(Region_number,Population,Variable) %>%
  summarize(avg=mean(Value),N=n()) %>%
  ungroup()

bad_regions <- pattern_summary %>%
  group_by(Region_number) %>%
  summarize(N=mean(N)) %>%
  ungroup() %>%
  filter(N <= 5)


good_regions <- pattern_summary %>% 
  spread(key= "Variable", value="avg") %>%
  mutate(out = Fst >= 1 & (pi >= 0.1 | dxy >= 1),out = ifelse(is.na(out),FALSE,out)) %>%
  group_by(Region_number) %>%
  summarize(not_match = sum(out)) %>%
  ungroup() %>%
  filter(not_match == 0)


interesting_speciation <- pattern_summary %>% 
  spread(key= "Variable", value="avg") %>%
  mutate(out = Fst >= 1 & dxy >= 1,out = ifelse(is.na(out),FALSE,out)) %>%
  group_by(Region_number) %>%
  summarize(not_match = sum(out)) %>%
  ungroup() %>%
  filter(not_match > 0,
         !(Region_number %in% filter_df$Region_number))

speciation_coastal <- speciation_df_single %>%
  filter(Region_number %in% interesting_speciation$Region_number)

ggplot(filter(pattern_summary,Region_number %in% interesting_speciation$Region_number),
       aes(x=Population,y=avg,colour=Variable))+
  geom_point()+
  facet_wrap(~Region_number)+
  theme(axis.text.x=element_text(angle=90))

speciation_df_single <- speciation_df_single %>% 
  filter(!(Region_number %in% bad_regions$Region_number),
         Region_number %in% good_regions$Region_number, 
         !(Region_number %in% filter_df$Region_number))

### Unifying neighbouring regions
speciation_df_single$Region_number <- 0
x <- 1
for (i in 1:(nrow(speciation_df_single)-1)){
  speciation_df_single$Region_number[i] <- x
  if(speciation_df_single$end[i]!=(speciation_df_single$start[i+1]-1)){
    x <- x+1
  }
}
speciation_df_single$Region_number[nrow(speciation_df_single)] <- speciation_df_single$Region_number[nrow(speciation_df_single)-1]+1 

speciation_df <- speciation_df_single %>%
  group_by(Chromosome,Region_number) %>%
  summarize(start = min(start),end=max(end),mid=mean(mid),sites=sum(sites),
            region_start = ((start+end)/2)-250000,
            region_end = ((start+end)/2)+250000,
            pi_LA1963=mean(pi_LA1963),pi_LA2931=mean(pi_LA2931),pi_LA2932=mean(pi_LA2932),
            pi_LA3111=mean(pi_LA3111),pi_LA4107=mean(pi_LA4107),pi_LA4330=mean(pi_LA4330),
            pi_LA_Peruvianum=mean(pi_LA_Peruvianum),
            dxy_LA1963_LA2931=mean(dxy_LA1963_LA2931),
            dxy_LA1963_LA2932=mean(dxy_LA1963_LA2932),
            dxy_LA1963_LA3111=mean(dxy_LA1963_LA3111),
            dxy_LA1963_LA4107=mean(dxy_LA1963_LA4107),
            dxy_LA1963_LA4330=mean(dxy_LA1963_LA4330),
            dxy_LA1963_LA_Peruvianum=mean(dxy_LA1963_LA_Peruvianum),
            dxy_LA1963_PI_Peruvianum=mean(dxy_LA1963_PI_Peruvianum),
            Dxy=mean(Dxy),Fst=mean(Fst), ZFst=mean(ZFst),ZDxy=mean(ZDxy),
            Fst_LA1963_LA2931=mean(Fst_LA1963_LA2931),
            Fst_LA1963_LA2932=mean(Fst_LA1963_LA2932),
            Fst_LA1963_LA3111=mean(Fst_LA1963_LA3111),
            Fst_LA1963_LA4107=mean(Fst_LA1963_LA4107),
            Fst_LA1963_LA4330=mean(Fst_LA1963_LA4330),
            Fst_LA1963_LA_Peruvianum=mean(Fst_LA1963_LA_Peruvianum),
            Fst_LA1963_PI_Peruvianum=mean(Fst_LA1963_PI_Peruvianum),
            Z_pi_LA1963=mean(Z_pi_LA1963),Z_pi_LA_Peruvianum=mean(Z_pi_LA_Peruvianum)) %>%
  ungroup() %>%
  mutate(window_size=end-start)

speciation_df %>%
  mutate(start = start/1000000,end=end/10000000) %>%
  select(Chromosome,start,end) %>%
  xtable(floating=FALSE,latex.environments=NULL,
         booktabs=TRUE,label="Tab:S_windows")

speciation_df %>%
  group_by(Chromosome) %>%
  summarize(Coverage=sum(window_size+1)/1000) %>%
  ggplot(aes(Chromosome,Coverage))+
  geom_col(fill="black")+
  theme_minimal()+
  ylab("Coverage in kb")+
  theme(axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

fwrite(speciation_df, file="islands_of_divergence.csv")
