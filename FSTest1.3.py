print('**************************************************************')
print('*                 FSTest software v 1.3                      *')
print('*       Seyed Milad Vahedi, Siavash Salek Ardestani          *')
print('*         smvahedi@dal.ca, siasia6650@gmail.com              *')
print('**************************************************************')
import time
import warnings
warnings.filterwarnings("ignore")
import allel
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import savefig
import io
import os

plt.switch_backend('agg')

FST = argparse.ArgumentParser()
FST.add_argument('--pop1',help='Input VCF file of population 1', required=True)
FST.add_argument('--pop2',help='Input VCF file of population 2', required=True)
FST.add_argument('--m',help='Fst estimation method: 1.Hudson , 2.Nei, 3.Weir&Cockerham, 4.Wright', type=int, required=True)
FST.add_argument('--zt',help='Fst Z-transformation: 1.SNP-based, 2.Win-based (optional)', type=int, required=False,default=0)
FST.add_argument('--win',help='Window size (optional)', type=int,required=False,default=0)
FST.add_argument('--step',help='Step size (optional)', type=int,required=False,default=0)
FST.add_argument('--mp',help='Fst Manhattan plot: 1.SNP-based, 2.Win-based (optional)', type=int,required=False,default=0)
FST.add_argument('--ztmp',help='Manhattan plot of Z-transformed Fst: 1.SNP-based, 2.Win-based (optional)', type=int,required=False,default=0)
FST.add_argument('--sl',help='Manhattan plot suggestive line (optional)', type=float,required=False,default=0.01)
FST.add_argument('--dpi',help='Plot dpi (optional)', type=float,required=False,default=300)
FST.add_argument('--o',help='Output files prefix', type=str, required=True)
args = FST.parse_args()

start_snp = time.time()

pop1 = allel.read_vcf(args.pop1)
pop2 = allel.read_vcf(args.pop2)

gt1 = allel.GenotypeArray(pop1['calldata/GT'])
gt2 = allel.GenotypeArray(pop2['calldata/GT'])

ac1 = gt1.count_alleles()
ac1 = pd.DataFrame(ac1)
ac1.columns=["np1","nq1"]
n1=len(pop1['samples'])

ac2 = gt2.count_alleles()
ac2 = pd.DataFrame(ac2)
ac2.columns=["np2","nq2"]
n2=len(pop2['samples'])

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    
snp = read_vcf(args.pop1)   
snp=snp[['CHROM','POS','ID']]

map = pd.concat([snp, ac1, ac2], axis=1)

map['AC1'] = map['np1'] + map['nq1']
map['AC2'] = map['np2'] + map['nq2']
map['p1'] = map['np1'] / map['AC1']
map['q1'] = 1-map['p1']
map['p2'] = map['np2'] / map['AC2']
map['q2'] = 1-map['p2']
map['p'] = (map['np1'] + map['np2'])/( map['AC2']+ map['AC1'])
map['q'] = 1- map['p']

if args.m == 1:
  X="Hud"
  map['Hud']=(((map['p1']-map['p2'])**2)-((map['p1']*map['q1'])/(n1-1))-((map['p2']*map['q2'])/(n2-1)))/((map['p1']*map['q2'])+(map['p2']*map['q1']))

if args.m == 2:
  X="Nei"
  map['Nei']=((map['p1']-map['p2'])**2)/((((map['p1']+map['p2'])/2)*2)*(1-((map['p1']+map['p2'])/2)))

if args.m == 3:
  X="WC"
  map['WC']=1-((2*((n1*n2)/(n1+n2))*(1/(n1+n2-2))*(((n1*map['p1'])*map['q1'])+(n2*map['p2']*map['q2'])))/((((n1*n2)/(n1+n2))*(((map['p1']-map['p2'])**2)))+(((2*((n1*n2)/(n1+n2)))-1)*(1/(n1+n2-2))*(((n1*map['p1'])*map['q1'])+(n2*map['p2']*map['q2'])))))                                        

if args.m == 4:
  X="Wright"
  map['Wright']=((2*map['p']*map['q'])-(((n1/(n1+n2))*2*map['p1']*map['q1'])+((n2/(n1+n2))*2*map['p2']*map['q2'])))/(2*map['p']*map['q'])

  
map1=map[['CHROM','POS', 'ID', X]]



if args.zt==1:
  map['Z(Fst)']=(map[X]-map[X].mean())/map[X].std()
  map3=map[['CHROM','POS', 'ID', 'Z(Fst)']]
  map1= map1.merge(map3, on=['CHROM','POS', 'ID'])

if args.zt==2:
  map['Z(Fst)']=(map[X]-map[X].mean())/map[X].std()
  map3=map[['CHROM','POS', 'ID', 'Z(Fst)']]
  map1= map1.merge(map3, on=['CHROM','POS', 'ID'])
  

Out=map1.reset_index(drop=True)
Out.to_csv(args.o+'.snp', sep='\t',index=False)

end_snp = time.time()

print("Elapsed time for processing SNPs:","%.2f" %round (end_snp-start_snp, 2), "s")


if args.mp==1:
    nsnp=len(Out)
    sugg1=int(args.sl*nsnp)
    sugg2=Out.sort_values(([X]), ascending=[False])
    sugg3=sugg2[X].iloc[sugg1]
    sugg4="{:.2f}".format(sugg3)
    print("Suggestive line of Manhattan plot of SNP-based Fst:",sugg4)
    Out['ind'] = range(len(Out))
    Out_grouped = Out.groupby(('CHROM'))
    fig_fstsnp = plt.figure()
    ax= plt.gca()
    colors = ['gray','black']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(Out_grouped):
        group.plot(kind='scatter', x='ind', y=X ,color=colors[num % len(colors)],  ax=ax,s= 0.3)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(Out)])
    ax.set_ylim([0,(Out[X].max()+0.2)])
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('FST ('+X+')')
    ax.tick_params(axis='both', which='major', labelsize=4)
    plt.axhline(y=sugg3, color='red', linestyle='--', linewidth = 0.9)
    fig_fstsnp.savefig(args.o+'.snpplot.png', dpi=args.dpi,format="png")
    
    
    
if args.ztmp==1:
    nsnp=len(Out)
    sugg1=int(args.sl*nsnp)
    sugg2=Out.sort_values((['Z(Fst)']), ascending=[False])
    sugg3=sugg2['Z(Fst)'].iloc[sugg1]
    sugg4="{:.2f}".format(sugg3)
    print("Suggestive line of Manhattan plot of SNP-based Z(Fst):",sugg4)
    Out['ind'] = range(len(Out))
    Out_grouped = Out.groupby(('CHROM'))
    fig_fstsnp = plt.figure()
    ax= plt.gca()
    colors = ['gray','black']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(Out_grouped):
        group.plot(kind='scatter', x='ind', y='Z(Fst)' ,color=colors[num % len(colors)],  ax=ax,s= 0.3)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(Out)])
    ax.set_ylim([0,(Out['Z(Fst)'].max()+0.2)])
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('ZFst('+X+')')
    ax.tick_params(axis='both', which='major', labelsize=4)
    plt.axhline(y=sugg3, color='red', linestyle='--', linewidth = 0.9)
    fig_fstsnp.savefig(args.o+'.snpplot.Z(Fst).png', dpi=args.dpi,format="png")


window = args.win
step = args.step
if window==0 or step == 0 :print >> exit()
start_win = time.time()

import pandas as pd
map2=pd.read_csv(args.o+'.snp',delimiter=r"\s+",header=0)

if args.m == 1:
  X="Hud"
if args.m == 2:
  X="Nei"
if args.m == 3:
  X="WC"
if args.m == 4:
  X="Wright"

df=map2[['CHROM','POS',X]]
indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=window)
df2 = df.join(df.POS.shift(-(window-1)), lsuffix='_start', rsuffix='_end')
df2 = df2.assign(Mean=df2.pop(X).rolling(window=indexer).mean()).iloc[::step]
df2 = df2[df2.POS_start.lt(df2.POS_end)].dropna()
df2['POS_end'] = df2['POS_end'].astype(int)
df2.columns = ['CHROM', 'START','END',X]


if args.zt==2:
    df3=map2[['CHROM','POS', 'Z(Fst)']]
    indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=window)
    df4 = df3.join(df3.POS.shift(-(window-1)), lsuffix='_start', rsuffix='_end')
    df4 = df4.assign(Mean=df4.pop('Z(Fst)').rolling(window=indexer).mean()).iloc[::step]
    df4 = df4[df4.POS_start.lt(df4.POS_end)].dropna()
    df4['POS_end'] = df4['POS_end'].astype(int)
    df4.columns = ['CHROM', 'START','END','Z(Fst)']
    df2 = df4.merge(df2, on=['CHROM', 'START','END'])

df2.to_csv(args.o+'.win', sep='\t',index=False)

end_win = time.time()

print("Elapsed time for processing windows:","%.2f" %round (end_win-start_win, 2), "s")


if args.mp==2:
    nsnp=len(df2)
    sugg1=int(args.sl*nsnp)
    sugg2=df2.sort_values(([X]), ascending=[False])
    sugg3=sugg2[X].iloc[sugg1]
    sugg4="{:.2f}".format(sugg3)
    print("Suggestive line of Manhattan plot of win-based Fst:",sugg4)
    df2['ind'] = range(len(df2))
    df2_grouped = df2.groupby(('Chr'))
    fig_fstsnp = plt.figure()
    ax= plt.gca()
    colors = ['gray','black']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df2_grouped):
        group.plot(kind='scatter', x='ind', y=X ,color=colors[num % len(colors)],  ax=ax,s= 0.3)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(df2)])
    ax.set_ylim([0,(df2[X].max()+0.2)])
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('FST ('+X+')')
    ax.tick_params(axis='both', which='major', labelsize=4)
    plt.axhline(y=sugg3, color='red', linestyle='--', linewidth = 0.9)
    fig_fstsnp.savefig(args.o+'.winplot.png', dpi=args.dpi,format="png")

    
if args.ztmp==2:
    nsnp=len(df2)
    sugg1=int(args.sl*nsnp)
    sugg2=df2.sort_values((['Z(Fst)']), ascending=[False])
    sugg4=sugg2['Z(Fst)'].iloc[sugg1]
    sugg5="{:.2f}".format(sugg4)
    print("Suggestive line of Manhattan plot of win-based Z(Fst):",sugg5)
    df2['ind'] = range(len(df2))
    df2_grouped = df2.groupby(('Chr'))
    fig_fstsnp = plt.figure()
    ax= plt.gca()
    colors = ['gray','black']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df2_grouped):
        group.plot(kind='scatter', x='ind', y='Z(Fst)' ,color=colors[num % len(colors)],  ax=ax,s= 0.3)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(df2)])
    ax.set_ylim([0,(df2['Z(Fst)'].max()+0.2)])
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('ZFst('+X+')')
    ax.tick_params(axis='both', which='major', labelsize=4)
    plt.axhline(y=sugg4, color='red', linestyle='--', linewidth = 0.9)
    fig_fstsnp.savefig(args.o+'.winplot.Z(Fst).png', dpi=args.dpi,format="png")



