# FSTest
<b>Estimation of different cross-population fixation indexes using parallel computing</b>
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
