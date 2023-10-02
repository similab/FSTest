![Logo](https://user-images.githubusercontent.com/52033302/116609032-435f0580-a90a-11eb-8d19-b791e713918d.png)

<b>Estimation of multiple cross-population fixation indexes using parallel computing</b>
</br>

<b> Description</b>
</br>
FSTest is a Python script for identifying and visualizing SNPs, or genomic windows loci contributing to the population differentiation using the Fst statistics described by Hudson, Nei, Weir & Cockerham, and Wright. The software supports the Z-transformation introduced by Akey. FSTest can handle missing genotypes and works with both phased and unphased data. The Fst statistics can be estimated per SNPs, or windows that are determined by a fixed number of SNPs. 
</br>


</pre></br>
<b> Parameter description</b>
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
