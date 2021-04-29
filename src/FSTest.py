print('**************************************************************')
print('*                 FSTest software v 1.0                      *')
print('*  Siavash Salek Ardestani, Seyed Milad Vahedi, Younes Miar  *')
print('*                  siasia6650@gmail.com                      *')
print('**************************************************************')

import warnings
warnings.filterwarnings("ignore")
import modin.pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import savefig
plt.switch_backend('agg')


FST = argparse.ArgumentParser()
FST.add_argument('--vcf',help='Input SNP VCF Format', required=True)
FST.add_argument('--g',help='Root file name of the groups file', required=True)
FST.add_argument('--chr',help='Number of chromosomes in the VCF file', type=int, required=True)
FST.add_argument('--n',help='Number of samples in the VCF file', type=int, required=True)
FST.add_argument('--m',help='FST estimation method: 1.Hudson , 2.Nei, 3.Weir&Cockerham, 4.Wright', type=int, required=True)
FST.add_argument('--di',help='di estiamtion of FST (Akey): 1.SNP-based, 2.Win-based (optional)', type=int, required=False,default=0)
FST.add_argument('--win',help='Window size (optional)', type=int,required=False,default=0)
FST.add_argument('--step',help='Step size (optional)', type=int,required=False,default=0)
FST.add_argument('--mp',help='FST Manhattan plot: 1.SNP-based, 2.Win-based (optional)', type=int,required=False,default=0)
FST.add_argument('--dimp',help='di Manhattan plot: 1.SNP-based, 2.Win-based (optional)', type=int,required=False,default=0)
FST.add_argument('--sl',help='Manhattan plot suggestive line (optional)', type=float,required=False,default=0.01)
FST.add_argument('--dpi',help='Plot dpi (optional)', type=float,required=False,default=300)
FST.add_argument('--o',help='Output files prefix', type=str, required=True)
args = FST.parse_args()

Chrs=args.chr
skip=Chrs+5
nind=args.n
shape=nind+9
VCF=pd.read_csv(args.vcf,delimiter=r"\s+",header=0,skiprows=skip,usecols=range(9,shape))
map=pd.read_csv(args.vcf,delimiter=r"\s+",header=0,skiprows=skip,usecols=range(0,3))
ID_Group=pd.read_csv(args.g,delimiter=r"\s+",header=0,dtype={'ID': str, 'Group': int})
map.columns=["Chr","Position","SNP"]
group1= ID_Group.loc[ID_Group.Group==1]
lst1 = group1['ID']
group2= ID_Group.loc[ID_Group.Group==2]
lst2 = group2['ID']
n1=len(lst1)
n2=len(lst2)
snp1=VCF[VCF.columns.intersection(lst1)]
snp2=VCF[VCF.columns.intersection(lst2)]
map[['np1', 'nq1']] = snp1.stack().str.extractall(r'(0)|(1)').notna().sum(level=0)
map[['np2', 'nq2']] = snp2.stack().str.extractall(r'(0)|(1)').notna().sum(level=0)
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

  
map1=map[['Chr','Position', 'SNP', X]]



if args.di==1:
  map['di']=(map[X]-map[X].mean())/map[X].std()
  map3=map[['Chr','Position', 'SNP', 'di']]
  map1= map1.merge(map3, on=['Chr','Position', 'SNP'])

if args.di==2:
  map['di']=(map[X]-map[X].mean())/map[X].std()
  map3=map[['Chr','Position', 'SNP', 'di']]
  map1= map1.merge(map3, on=['Chr','Position', 'SNP'])
  

Out=map1.reset_index(drop=True)
Out.to_csv(args.o+'.snp', sep='\t',index=False)

if args.mp==1:
    nsnp=len(Out)
    sugg1=int(args.sl*nsnp)
    sugg2=Out.sort_values(([X]), ascending=[False])
    sugg3=sugg2[X].iloc[sugg1]
    print("Suggestive line of SNP-based FST Manhattan plot:",sugg3)
    Out['ind'] = range(len(Out))
    Out_grouped = Out.groupby(('Chr'))
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
    
    
    
if args.dimp==1:
    nsnp=len(Out)
    sugg1=int(args.sl*nsnp)
    sugg2=Out.sort_values((['di']), ascending=[False])
    sugg3=sugg2['di'].iloc[sugg1]
    print("Suggestive line of SNP-based di Manhattan plot:",sugg3)
    Out['ind'] = range(len(Out))
    Out_grouped = Out.groupby(('Chr'))
    fig_fstsnp = plt.figure()
    ax= plt.gca()
    colors = ['gray','black']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(Out_grouped):
        group.plot(kind='scatter', x='ind', y='di' ,color=colors[num % len(colors)],  ax=ax,s= 0.3)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(Out)])
    ax.set_ylim([0,(Out['di'].max()+0.2)])
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('di ('+X+')')
    ax.tick_params(axis='both', which='major', labelsize=4)
    plt.axhline(y=sugg3, color='red', linestyle='--', linewidth = 0.9)
    fig_fstsnp.savefig(args.o+'.snpplot.di.png', dpi=args.dpi,format="png")


window = args.win
step = args.step
if window==0 or step == 0 :print >> exit()
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

df=map2[['Chr','Position',X]]
indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=window)
df2 = df.join(df.Position.shift(-(window-1)), lsuffix='_start', rsuffix='_end')
df2 = df2.assign(Mean=df2.pop(X).rolling(window=indexer).mean()).iloc[::step]
df2 = df2[df2.Position_start.lt(df2.Position_end)].dropna()
df2['Position_end'] = df2['Position_end'].astype(int)
df2.columns = ['Chr', 'Start_pos','End_pos',X]


if args.di==2:
    df3=map2[['Chr','Position', 'di']]
    indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=window)
    df4 = df3.join(df3.Position.shift(-(window-1)), lsuffix='_start', rsuffix='_end')
    df4 = df4.assign(Mean=df4.pop('di').rolling(window=indexer).mean()).iloc[::step]
    df4 = df4[df4.Position_start.lt(df4.Position_end)].dropna()
    df4['Position_end'] = df4['Position_end'].astype(int)
    df4.columns = ['Chr', 'Start_pos','End_pos','di']
    df2 = df4.merge(df2, on=['Chr', 'Start_pos','End_pos'])

df2.to_csv(args.o+'.win', sep='\t',index=False)


if args.mp==2:
    nsnp=len(df2)
    sugg1=int(args.sl*nsnp)
    sugg2=df2.sort_values(([X]), ascending=[False])
    sugg3=sugg2[X].iloc[sugg1]
    print("Suggestive line of win-based FST Manhattan plot:",sugg3)
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

    
if args.dimp==2:
    nsnp=len(df2)
    sugg1=int(args.sl*nsnp)
    sugg2=df2.sort_values((['di']), ascending=[False])
    sugg4=sugg2['di'].iloc[sugg1]
    print("Suggestive line of win-based di Manhattan plot:",sugg4)
    df2['ind'] = range(len(df2))
    df2_grouped = df2.groupby(('Chr'))
    fig_fstsnp = plt.figure()
    ax= plt.gca()
    colors = ['gray','black']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df2_grouped):
        group.plot(kind='scatter', x='ind', y='di' ,color=colors[num % len(colors)],  ax=ax,s= 0.3)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(df2)])
    ax.set_ylim([0,(df2['di'].max()+0.2)])
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('di ('+X+')')
    ax.tick_params(axis='both', which='major', labelsize=4)
    plt.axhline(y=sugg4, color='red', linestyle='--', linewidth = 0.9)
    fig_fstsnp.savefig(args.o+'.winplot.di.png', dpi=args.dpi,format="png")
    