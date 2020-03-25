import pandas as pd
import numpy as np
import argparse as parser
pd.options.mode.chained_assignment = None

# Reading in the arguments
parser = parser.ArgumentParser(usage= 'This program searches for genomic regions with a higher differentiation than expected under neutrality between two diverged species. If other populations of one species are given a second filter can be applied to see that the regions do not differ within one species.')

parser.add_argument('--input',action='store',help='Input file created with popgenWindows.py found at: https://github.com/simonhmartin/genomics_general',required= True,metavar='Input')
parser.add_argument('--pop1',action='store',help='Name of the first population, which is genetically closest to the second population',required= True,metavar='pop1')
parser.add_argument('--pop2',action='store',help='Name of the second population from the species to which the differentiation is measured',required= True,metavar='pop2')
parser.add_argument('--pop3',action='store',help='Name of the third population, which is genetically distant to the second population and the same species as the first population',required= False,metavar='pop3')
parser.add_argument('--pop3_list',action='store',help='Same as "pop3", but takes a file with more populations. One population has to be given per line',required= False,metavar='pop3_list')
parser.add_argument('--pop4',action='store',help='Name of the fourth population, which is genetically distant to the first population and the same species as the second population',required= False,metavar='pop4')
parser.add_argument('--pop4_list',action='store',help='Same as "pop4", but takes a file with more populations. One population has to be given per line',required= False,metavar='pop4_list')
parser.add_argument('--output',action='store',help='Output prefix. Default = output',required= False, metavar='Output',default='output')
parser.add_argument('--thres',action='store',help='Threshold for the outlier scan. Defines how many standard deviations above the mean the dxy and Fst should be. Default = 2',required= False,metavar='Threshold', default=2,type=int)
parser.add_argument('--min_snps',action='store',help='Minimum amount of SNPs per window. Windows below the count will be filtered out. Default = 20',required= False,metavar='Minimum_SNPs_per_window', default=20,type=int)


args = parser.parse_args()

pop1 = args.pop1
pop2 = args.pop2

if args.pop3 is not None:
    pop3_list = [args.pop3]
elif args.pop3_list is not None:
    with open(args.pop3_list) as f:
        pop3_list = f.read().splitlines()
else:
    pop3_list = None

if args.pop4 is not None:
    pop4_list = [args.pop4]
elif args.pop4_list is not None:
    with open(args.pop4_list) as f:
        pop4_list = f.read().splitlines()
else:
    pop4_list = None  
    
popstat = pd.read_csv(args.input)
n_chr = len(popstat.scaffold.unique())
chr_new = []
for i in range(0,n_chr):
    chr_new.append('Chr' + str(i+1))

chr_old = popstat.scaffold.unique()
chr_dict = dict(zip(chr_old,chr_new))
popstat['scaffold'] = popstat.scaffold.map(chr_dict)
popstat['scaffold'] = popstat.scaffold.astype('category')
popstat = popstat[popstat.sites >= args.min_snps]
    
# Start run with only pop1 and pop2    
if ((pop3_list is None) & (pop4_list is None)):    
    popstat = pd.read_csv(args.input)
    n_chr = len(popstat.scaffold.unique())
    chr_new = []
    for i in range(0,n_chr):
        chr_new.append('Chr' + str(i+1))

    chr_old = popstat.scaffold.unique()
    chr_dict = dict(zip(chr_old,chr_new))
    popstat['scaffold'] = popstat.scaffold.map(chr_dict)
    popstat['scaffold'] = popstat.scaffold.astype('category')
    popstat = popstat[popstat.sites >= args.min_snps]

# Filter out the relevant columns and z-transform them
    pre_filter = popstat.iloc[:,(popstat.columns.str.contains(pop1) & np.logical_not(popstat.columns.str.contains('pi')) & popstat.columns.str.contains(pop2))]

    pre_filter = pre_filter.apply(lambda x: (x - x.mean())/x.std(ddof=0))

# Preparing the filter for the outlier
    filter_inter = pre_filter >= args.thres
    filter_inter['Keep'] = filter_inter.apply(lambda x: sum(x) == 2, axis = 1)
    filter_list = filter_inter.Keep
    islands_single = popstat.iloc[filter_list.to_list(),:]

# Start with pop1, pop2, and pop3
elif ((pop4_list is None) & (pop3_list is not None)):
    pre_filter = popstat.iloc[:,(popstat.columns.str.contains(pop1) & np.logical_not(popstat.columns.str.contains('pi')))]
    
    pop3_columns = pd.DataFrame(columns=pop3_list)
    for i in pop3_list:
        pop3_columns[i] = pre_filter.columns.str.contains(i)
        
    pop3_columns['Keep'] = pop3_columns.sum(axis=1)
    pop3_columns['Keep'] = pop3_columns.Keep + pre_filter.columns.str.contains(pop2)
    pop3_columns['Keep'] = pop3_columns.Keep > 0
    pre_filter = pre_filter.iloc[:,pop3_columns.Keep.to_list()]
    
    pre_filter = pre_filter.apply(lambda x: (x - x.mean())/x.std(ddof=0))
    
    filter_inter = pre_filter.iloc[:,pre_filter.columns.str.contains(pop2)] >= args.thres
    filter_inter['Keep'] = filter_inter.apply(lambda x: sum(x) == 2, axis = 1)
    filter_intra = pre_filter.iloc[:,np.logical_not(pre_filter.columns.str.contains(pop2))&pre_filter.columns.str.contains('dxy')] >= 0
    filter_intra['n_diff'] = filter_intra.sum(axis=1)
    filter_intra['Keep'] = filter_intra.n_diff == 0
    filter_list = pd.concat([filter_intra.Keep,filter_inter.Keep],axis=1)
    filter_list['outlier'] = filter_list.apply(lambda x: sum(x) == 2,axis =1)
    islands_single = popstat.iloc[filter_list.outlier.to_list(),:]

# Start with pop1,pop2, and pop4 here
elif ((pop3_list is None) & (pop4_list is not None)):
    pre_filter = popstat.iloc[:,(popstat.columns.str.contains(pop2) & np.logical_not(popstat.columns.str.contains('pi')))]
    pop4_columns = pd.DataFrame(columns=pop4_list)
    for i in pop4_list:
        pop4_columns[i] = pre_filter.columns.str.contains(i)
        
    pop4_columns['Keep'] = pop4_columns.sum(axis=1)
    pop4_columns['Keep'] = pop4_columns.Keep + pre_filter.columns.str.contains(pop1)
    pop4_columns['Keep'] = pop4_columns.Keep > 0
    pre_filter = pre_filter.iloc[:,pop4_columns.Keep.to_list()]
    
    pre_filter = pre_filter.apply(lambda x: (x - x.mean())/x.std(ddof=0))
    
    filter_inter = pre_filter.iloc[:,pre_filter.columns.str.contains(pop1)] >= args.thres
    filter_inter['Keep'] = filter_inter.apply(lambda x: sum(x) == 2, axis = 1)
    filter_intra = pre_filter.iloc[:,np.logical_not(pre_filter.columns.str.contains(pop1))&pre_filter.columns.str.contains('dxy')] >= 0
    filter_intra['n_diff'] = filter_intra.sum(axis=1)
    filter_intra['Keep'] = filter_intra.n_diff == 0
    filter_list = pd.concat([filter_intra.Keep,filter_inter.Keep],axis=1)
    filter_list['outlier'] = filter_list.apply(lambda x: sum(x) == 2,axis =1)
    islands_single = popstat.iloc[filter_list.outlier.to_list(),:]

# Start with pop1,pop2,pop3, and pop4 here
elif ((pop4_list is not None) & (pop3_list is not None)):

## Preparation of pop1 filtering
    pre_filter_pop1 = popstat.iloc[:,(popstat.columns.str.contains(pop1) & np.logical_not(popstat.columns.str.contains('pi')))]
    pop3_columns = pd.DataFrame(columns=pop3_list)
    for i in pop3_list:
        pop3_columns[i] = pre_filter_pop1.columns.str.contains(i)
    pop3_columns['Keep'] = pop3_columns.sum(axis=1)
    pop3_columns['Keep'] = pop3_columns.Keep + pre_filter_pop1.columns.str.contains(pop2)
    pop3_columns['Keep'] = pop3_columns.Keep > 0
    pre_filter_pop1 = pre_filter_pop1.iloc[:,pop3_columns.Keep.to_list()]    
    pre_filter_pop1 = pre_filter_pop1.apply(lambda x: (x - x.mean())/x.std(ddof=0))

## Preparation of pop2 filtering
    pre_filter_pop2 = popstat.iloc[:,(popstat.columns.str.contains(pop2) & np.logical_not(popstat.columns.str.contains('pi')))]
    pop4_columns = pd.DataFrame(columns=pop4_list)
    for i in pop4_list:
        pop4_columns[i] = pre_filter_pop2.columns.str.contains(i)
        
    pop4_columns['Keep'] = pop4_columns.sum(axis=1)
    pop4_columns['Keep'] = pop4_columns.Keep + pre_filter_pop2.columns.str.contains(pop1)
    pop4_columns['Keep'] = pop4_columns.Keep > 0
    pre_filter_pop2 = pre_filter_pop2.iloc[:,pop4_columns.Keep.to_list()]
    
    pre_filter_pop2 = pre_filter_pop2.apply(lambda x: (x - x.mean())/x.std(ddof=0))
    
## Putting filter together
    filter_inter = pre_filter_pop1.iloc[:,pre_filter_pop1.columns.str.contains(pop2)] >= args.thres
    filter_inter['Keep'] = filter_inter.apply(lambda x: sum(x) == 2, axis = 1)
    
    filter_intra_pop1 = pre_filter_pop1.iloc[:,np.logical_not(pre_filter_pop1.columns.str.contains(pop2))&pre_filter_pop1.columns.str.contains('dxy')] >= 0
    filter_intra_pop1['n_diff'] = filter_intra_pop1.sum(axis=1)
    filter_intra_pop1['Keep'] = filter_intra_pop1.n_diff == 0
 
    filter_intra_pop2 = pre_filter_pop2.iloc[:,np.logical_not(pre_filter_pop2.columns.str.contains(pop1))&pre_filter_pop2.columns.str.contains('dxy')] >= 0
    filter_intra_pop2['n_diff'] = filter_intra_pop2.sum(axis=1)
    filter_intra_pop2['Keep'] = filter_intra_pop2.n_diff == 0

    filter_list = pd.concat([filter_intra_pop1.Keep,filter_intra_pop2.Keep,filter_inter.Keep],axis=1)
    filter_list['outlier'] = filter_list.apply(lambda x: sum(x) == 3,axis =1)
    islands_single = popstat.iloc[filter_list.outlier.to_list(),:]    

# Filtering out the single outlier windows and unifying neighbouring windows to islands
islands_single['Region_number'] = 0
x = 1
for i in range(0,(len(islands_single)-1)):
    islands_single.Region_number.iloc[i] = x
    if islands_single.iloc[i,:]['end'] != (islands_single.iloc[(i+1),:]['start'] -1):
        x = x+1

if (islands_single.end.iloc[-2] == (islands_single.start.iloc[-1]-1)):
    islands_single.Region_number.iloc[-1] = islands_single.Region_number.iloc[-2]
else:
    islands_single.Region_number.iloc[-1] = (islands_single.Region_number.iloc[-2]+1)

islands_single_sorted = islands_single.iloc[:,:5]
islands_single_sorted['Region_number'] = islands_single['Region_number']
    
islands_unified = islands_single_sorted.groupby(['scaffold','Region_number']).agg({'start':min,'end':max,'sites':sum}).dropna().sort_values('Region_number').reset_index()
islands_unified['window_size'] = islands_unified.end - islands_unified.start
islands_unified.columns.values[0]='Chromosome'
islands_unified.to_csv(args.output+'.csv',index=False)
