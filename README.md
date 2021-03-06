![Logo](https://user-images.githubusercontent.com/52033302/116609032-435f0580-a90a-11eb-8d19-b791e713918d.png)

<b>Estimation of multiple cross-population fixation indexes using parallel computing</b>
</br>
</br>
<b>1. Description</b>
</br>
FSTest is a Python script for identifying and visualizing genomic loci under selection using the "Fst" statistics described by Hudson, Nei, Weir & Cockerham, and Wright as well as Z-transformation introduced by Akey. The software employs parallel computing to reduce the calculation time and efforts have been made to add new features that would make it more beneficial for selection signatures analysis. FSTest can handle missing genotypes and works with both phased and unphased data.
</br>

<b>2. Download and installation</b>
</br>
<b>2.1. Download for linux</b>
</br>
<pre>
        git clone https://github.com/similab/FSTest.git;
	cd FSTest;
	cd src;
	chmod 775 FSTest.py
</pre>
</br>
<b>2.2. Dependencies installation</b>
 </br>
FSTest package works with python 3.8 or later versions and the following libraries are required:
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

```py
	Usage: python FSTest_v1.0.py --vcf  <in.vcf> --g <in.txt> --chr <int> --n <int> --m <int> --o <str>

               -h, --help   show this help message and exit
               --vcf VCF    Input SNP VCF Format
               --g G        Root file name of the groups file
               --chr CHR    Number of chromosomes in the VCF file
               --n N        Number of samples in the VCF file
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
<b>3. Examples</b>
</br>
Sample VCF and ID group input files are available in the "example" directory.</b>
 </br>
<b>3.1. SNP-based Fst estimation</b>
 </br>
```py
        python FSTest_v1.0.py --vcf sheep.vcf --g ID_Group.txt --chr 26 --n 133 --m 1 --zt 1 --mp 1 --ztmp 1 --sl 0.05 --dpi 600 --o test.snp
```
</br>
<b>Outputs</b>
</br>
1. Fst and Z(Fst) values of SNPs using selected method </br>
2. Manhattan plot of Z(Fst) values 
<img src="https://github.com/ymiarlab/FSTest/blob/main/result1/test.snp.snpplot.png" width="800"/>

<b>3.2. Window-based Fst estimation</b>
 </br>
```py
        python FSTest_v1.0.py --vcf sheep.vcf --g ID_Group.txt --chr 26 --n 133 --m 1 --zt 2 --win 20 --step 5 --mp 2 --ztmp 2 --sl 0.05 --dpi 600 --o test.win
```
<b>Outputs</b>
</br>
1. Fst and Z(Fst) values of windows using selected method </br>
2. Manhattan plot of Z(Fst) values
<img src="https://github.com/ymiarlab/FSTest/blob/main/result2/test.win.winplot.png" width="800"/>

<b>References</b>
1. Akey JM, Ruhe AL, Akey DT, Wong AK, Connelly CF, Madeoy J, Nicholas TJ, Neff MW. Tracking footprints of artificial selection in the dog genome. Proceedings of the National Academy of Sciences. 2010 Jan 19;107(3):1160-5.
2. Hudson RR, Slatkin M, Maddison WP. Estimation of levels of gene flow from DNA sequence data. Genetics. 1992 Oct 1;132(2):583-9.
3. Nei M. Definition and estimation of fixation indices. Evolution. 1986 May 1;40(3):643-5.
4. Weir BS, Cockerham CC. Estimating F-statistics for the analysis of population structure. Evolution. 1984 Nov 1:1358-70.
5. Wright S. The genetical structure of populations. Annals of Eugenics. 1949 Jan;15(1):323-54.

