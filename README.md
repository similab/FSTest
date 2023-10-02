![Logo](https://user-images.githubusercontent.com/52033302/116609032-435f0580-a90a-11eb-8d19-b791e713918d.png)

<b>Estimation of multiple cross-population fixation indexes using parallel computing</b>
</br>
<b>Citation:</b>
Salek Ardestani S, Vahedi SM. FSTest: an efficient tool for estimation and visualization of cross-population fixation index on variant call format files. bioRxiv. 2022 Jan 1;2022.08.16.503949. 
</br>

<b>1. Description</b>
</br>
FSTest is a Python script for identifying and visualizing genomic loci under selection using the "Fst" statistics described by Hudson, Nei, Weir & Cockerham, and Wright as well as Z-transformation introduced by Akey. The software employs parallel computing to reduce the calculation time and efforts have been made to add new features that would make it more beneficial for selection signatures analysis. FSTest can handle missing genotypes and works with both phased and unphased data.
</br>


</pre></br>
<b>2. Parameter description</b>
</br>

```py
	Usage: ./FSTest1.3 --pop1  <in.vcf> --pop2 <in.vcf> --m <int> --o <str>

               -h, --help   show this help message and exit
               --pop1 VCF   Input VCF file of population 1
               --pop2 VCF   Input VCF file of population 2
               --m M        FST estimation method: 1.Hudson , 2.Nei, 3.Weir&Cockerham, 4.Wright
               --zt ZT      Fst Z-transformation: 1.SNP-based, 2.Win-based (optional)
               --win WIN    Window size (optional)
               --step STEP  Step size (optional)
               --mp MP      FST Manhattan plot: 1.SNP-based, 2.Win-based (optional)
               --ztmp ZTMP  Manhattan plot of Z-transformed Fst: 1.SNP-based, 2.Win-based (optional)
               --sl SL      Manhattan plot suggestive line (optional)
               --dpi DPI    Plot dpi (optional)
               --o O        Output files prefix
```
