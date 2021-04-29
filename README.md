# FSTest
<b>Cross-population fixation index estimation using parallel computing</b>
</br>
<b>1. Description</b>
</br>
<b>2. Download and installation</b>
</br>
<pre>
        git clone https://github.com/Siavash-cloud/DEBV_calculator.git
        tar -zxvf  PopLDdecayXXX.tar.gz
</pre>
</br>
 <b>2.2. Dependencies installation</b>
 </br>
FSTest package works on python 3.8 or later versions.
</br>
<pre>
         pip install numpy
         pip install pandas
         pip install modin[all]
         pip install matplotlib
         pip install argparse
</pre></br>
<b>3. Parameter description</b>
</br>
```php
	Usage: python FSTest_v1.0.py --vcf  <in.vcf> --g <in.txt> --chr <int> --n <int> --m <int> --o <str>

               -h, --help   show this help message and exit
               --vcf VCF    Input SNP VCF Format
               --g G        Root file name of the groups file
               --chr CHR    Number of chromosomes in the VCF file
               --n N        Number of samples in the VCF file
               --m M        FST estimation method: 1.Hudson , 2.Nei, 3.Weir&Cockerham, 4.Wright
               --di DI      di estiamtion of FST (Akey): 1.SNP-based, 2.Win-based (optional)
               --win WIN    Window size (optional)
               --step STEP  Step size (optional)
               --mp MP      FST Manhattan plot: 1.SNP-based, 2.Win-based (optional)
               --dimp DIMP  di Manhattan plot: 1.SNP-based, 2.Win-based (optional)
               --sl SL      Manhattan plot suggestive line (optional)
               --dpi DPI    Plot dpi (optional)
               --o O        Output files prefix
```
</br>
Linkage disequilibrium (LD) decay[1] is the most important and most common analysis in the population resequencing[2]. Special in the self-pollinated crops, the LD decay may not only reveal much about domestication and breed history[3], but also can reveal gene flow phenomenon, selection regions[1].However, to measure the LD decay, it takes too much resources and time by using currently existent software and tools. The LD decay studies also generate extraordinarily large amounts of data to temporary storage when you using the mainstream software "Haploview"[4], the classical LD processing tools. Effective use and analysis to get the LD decay result remains a difficult task for individual researchers. Here, we introduce PopLDdecay, a simple- efficient software for LD decay analysis, which processes the Variant Call Format (VCF)[5] file to produce the LD decay statistics results and plot the LD decay graphs. PopLDdecay is designed to use compressed data files as input or output to save storage space and it facilitates faster and more computationally efficient than the currently existent softwares. This software makes the LD decay pipeline significantly
* <b> Parameter description</b>
```php
	Usage: PopLDdecay -InVCF  <in.vcf.gz>  -OutStat <out.stat>

		-InVCF       <str>    Input SNP VCF Format
		-InGenotype  <str>    Input SNP Genotype Format
		-OutStat     <str>    OutPut Stat Dist ~ r^2 File

		-SubPop      <str>    SubGroup SampleList of VCFFile [ALLsample]
		-MaxDist     <int>    Max Distance (kb) between two SNP [300]
		-MAF         <float>  Min minor allele frequency filter [0.005]
		-Het         <float>  Max ratio of het allele filter [0.88]
		-Miss        <float>  Max ratio of miss allele filter [0.25]
		-OutFilterSNP         OutPut the final SNP to calculate
		-OutPairLD   <int>    OutPut the PairWise SNP LD info [0]
		                      0/2:No_Out 1/3/4:Out_Brief 5:Out_Full

		-help                 Show more help [hewm2008 v3.30]

```

