![Logo](https://user-images.githubusercontent.com/52033302/116609032-435f0580-a90a-11eb-8d19-b791e713918d.png)

<b>Estimation of different cross-population fixation indexes using parallel computing</b>
</br>

Authors: Siavash Salek Ardestani, Seyed Milad Vahedi, Younes Miar
</br>

<b>1. Description</b>
</br>
FSTest is a script written in Python to identify and visualize genomic loci under selection using the "Fst" statistics described by Hudson, Nei, Weir & Cockerham, and Wright as well as "di" transformation introduced by Akey. The software employs parallel computing to boost the calculation process and efforts have being made to add new features that would make it more beneficial for selection signatures analysis.
</br>

<b>2. Download and installation</b>
</br>
<b>2.1. Download for linux</b>
</br>
<pre>
        git clone https://github.com/Miarlab/FSTest.git;
	cd FSTest;
	cd src;
	chmod 775 FSTest.py
</pre>
</br>
<b>2.2. Dependencies installation</b>
 </br>
FSTest package works on python 3.8 or later versions. Installing these dependencies prior to running FSTest are recommened:
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
               --di DI      di estiamtion of FST (Akey): 1.SNP-based, 2.Win-based (optional)
               --win WIN    Window size (optional)
               --step STEP  Step size (optional)
               --mp MP      FST Manhattan plot: 1.SNP-based, 2.Win-based (optional)
               --dimp DIMP  di Manhattan plot: 1.SNP-based, 2.Win-based (optional)
               --sl SL      Manhattan plot suggestive line (optional)
               --dpi DPI    Plot dpi (optional)
               --o O        Output files prefix
```
<b>3. Examples</b>
</br>
Sample VCF and ID group input files can be downloaded from the "example" directory.</b>
 </br>
<b>3.1. SNP-based Fst estimation</b>
 </br>
```py
        python FSTest_v1.0.py --vcf sheep.vcf --g ID_Group.txt --chr 26 --n 133 --m 1 --di 1 --mp 1 --dimp 1 --sl 0.05 --dpi 600 --o test.snp
```
</br>
<b>Outputs</b>
</br>
1. Fst and di values of SNPs using selected method (test.snp.snp)</br>
2. Manhattan plot of Fst values (snp.snpplot.di.png)
<img src="https://github.com/Miarlab/FSTest/blob/main/result1/test.snp.snpplot.png" width="800"/>



3. Manhattan plot of di values (snp.snpplot.di.png)
<img src="https://github.com/Miarlab/FSTest/blob/main/result1/test.snp.snpplot.di.png" width="800"/>

<b>3.2. Window-based Fst estimation</b>
 </br>
```py
        python FSTest_v1.0.py --vcf sheep.vcf --g ID_Group.txt --chr 26 --n 133 --m 1 --di 2 --win 20 --step 5 --mp 2 --dimp 2 --sl 0.05 --dpi 600 --o test.win
```
<b>Outputs</b>
</br>
1. Fst and di values of windows using selected method (test.win.win)</br>
2. Manhattan plot of Fst values (test.win.winplot.png)
<img src="https://github.com/Miarlab/FSTest/blob/main/result2/test.win.winplot.png" width="800"/>



3. Manhattan plot of di values (test.win.winplot.di.png)
<img src="https://github.com/Miarlab/FSTest/blob/main/result2/test.win.winplot.di.png" width="800"/>

<b>References</b>
1. Akey, Joshua M., et al. "Tracking footprints of artificial selection in the dog genome." Proceedings of the National Academy of Sciences 107.3 (2010): 1160-1165.
2. Hudson, Richard R., Montgomery Slatkin, and Wayne P. Maddison. "Estimation of levels of gene flow from DNA sequence data." Genetics 132.2 (1992): 583-589.
3. Nei, Masatoshi. "Definition and estimation of fixation indices." Evolution 40.3 (1986): 643-645.
4. Weir, Bruce S., and C. Clark Cockerham. "Estimating F-statistics for the analysis of population structure." evolution (1984): 1358-1370.
5. Wright, Sewall. "The genetical structure of populations." Annals of eugenics 15.1 (1949): 323-354.
