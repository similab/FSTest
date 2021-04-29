# FSTest
<b>Estimation of different cross-population fixation indexes using parallel computing</b>
</br>
<b>1. Description</b>
</br>
<b>2. Download and installation</b>
</br>
2.1. Download for linux</b>
</br>
<pre>
        git clone https://github.com/Miarlab/FSTest.git;
	cd FSTest;
	cd src;
	chmod 775 FSTest.py
</pre>
</br>
2.2. Dependencies installation</b>
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
Sample VCF and ID group input files can be downloaded form the "example" directory.</b>
 </br>
3.1. SNP-based Fst estimation</b>
 </br>
```py
        python FSTest_v1.0.py --vcf sheep.vcf --g ID_Group.txt --chr 26 --n 133 --m 1 --di 1 --mp 1 --dimp 1 --sl 0.05 --dpi 600 --o test.snp
```
</br>
<b>Outputs:</b>
</br>
1. Fst values of SNPs using selected method (test.snp.snp)</br>
2. Manhattan plot of Fst values</br>
![test snp snpplot](https://user-images.githubusercontent.com/52033302/116605840-51ab2280-a906-11eb-8f1f-f62d9d1d1547.png)



3. Manhattan plot of di values</br>
![test snp snpplot di](https://user-images.githubusercontent.com/52033302/116604963-32f85c00-a905-11eb-8839-0abccce2cbe2.png)

3.2. Window-based Fst estimation</b>
 </br>
```py
        python FSTest_v1.0.py --vcf sheep.vcf --g ID_Group.txt --chr 26 --n 133 --m 1 --di 2 --win 20 --step 5 --mp 2 --dimp 2 --sl 0.05 --dpi 600 --o test.win
```



