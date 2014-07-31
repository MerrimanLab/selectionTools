
<div class="item">

<h3>Documentation</h3>
<p>
The VCFtools package is broadly split into two sections:
<ul class="faq">
<li> The <span class="cmd">vcftools</span> binary program, generally used to analyse VCF files.
<li> The <span class="cmd">Vcf.pm</span> perl module, which is a general Perl API containing a core of
the utilities <span class="cmd">vcf-convert</span>, <span class="cmd">vcf-merge</span>, 
<span class="cmd">vcf-compare</span>, <span class="cmd">vcf-isec</span>, and others.
</ul>
</p>

<h3>Documentation</h3>
<ul class="faq">
<li> <a href="options.html">Documentation for the binary executable <span class="cmd">vcftools</span></a>
<li> <a href="perl_module.html">The Perl API and scripts</a>
<li> <a href="htslib.html">Fast htslib C tools</a> <span style="font-size:70%;font-style:italic;color:#ff5500;">(new)</span>
</ul>

<h3>Examples of usage by topic</h3>
<ul class="faq">
<li> <a href="#install">Installation</a>
<li> <a href="#annotate">Annotating</a>
<li> <a href="#compare">Comparing</a>
<li> <a href="#concat">Concatenating</a>
<li> <a href="#convert">Converting</a>
<li> <a href="#filter">Filtering</a>
<li> <a href="#isec">Intersections, complements</a>
<li> <a href="#merge">Merging</a>
<li> <a href="#query">Querying</a>
<li> <a href="#shuffle">Reordering columns</a>
<li> <a href="#stats">Stats</a>
<li> <a href="#subset">Stripping columns</a>
<li> <a href="#one-liners">Useful shell one-liners</a>
<li> <a href="#validator">Validating</a>
<li> <a href="#more">...and more</a>
</ul>
</div>


<div class="item">
<a name="install" class="Q">Installation</a>
<p> The VCFtools package can be decompressed by the command
<p class="codebox">tar -xzf vcftools_<i>version_number</i>_source.tar.gz</p>
<p> To build the <span class="cmd">vcftools</span> executable, type "make" in the vcftools folder. </p>
<p>
The <b>Perl scripts require that VCF files are compressed by <span class="cmd">bgzip</span>
and indexed by <span class="cmd">tabix</span></b> (both tools are part of the tabix package, 
available for <a href="https://sourceforge.net/projects/samtools/files/tabix/">download here</a>).
Both tools must be in directories that are listed in the <a href="http://en.wikipedia.org/wiki/PATH_%28variable%29">PATH</a> environment variable.
For running the Perl scripts, the <b>PERL5LIB environment variable must be set</b> to include the Vcf.pm
module</p>
<p class="codebox">export PERL5LIB=/path/to/your/vcftools-directory/perl</p>
The tools can be tested by running the script
<p class="codebox">/path/to/your/vcftools-directory/perl/test.t</p> 
If the command complains about missing <tt>Test::Most</tt> perl module, do not
worry, it is needed only for testing, not for running VCFtools.
</div>


<div class="item">
<a name="annotate" class="Q">Annotating</a>
<p></p>
<p class="codebox">
# Add custom annotations <br>
cat in.vcf | vcf-annotate -a annotations.gz \ <br>
&nbsp;&nbsp; -d key=INFO,ID=ANN,Number=1,Type=Integer,Description='My custom annotation' \ <br>
&nbsp;&nbsp; -c CHROM,FROM,TO,INFO/ANN &gt; out.vcf <br>
<br>
# Apply SnpCluster filter <br>
cat in.vcf | <a href="perl_module.html#vcf-annotate">vcf-annotate</a> --filter SnpCluster=3,10 &gt; out.vcf <br>
</p>
</div>



<div class="item">
<a name="compare" class="Q">Comparing</a>
<p></p>
<p class="codebox">
<a href="perl_module.html#vcf-compare">vcf-compare</a> A.vcf.gz B.vcf.gz C.vcf.gz <br>
<a href="htslib.html#check">vcf check</a> A.vcf.gz B.vcf.gz <br>
</p>
</div>


<div class="item">
<a name="concat" class="Q">Concatenating</a>
<p></p>
<p class="codebox"><a href="perl_module.html#vcf-concat">vcf-concat</a> A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz
</p>
</div>


<div class="item">
<a name="convert" class="Q">Converting</a>
<p></p>
<p class="codebox">
# Convert between VCF versions  <br>
zcat file.vcf.gz | <a href="perl_module.html#vcf-convert">vcf-convert</a> -r reference.fa | bgzip -c &gt; out.vcf.gz  <br>
<br>
# Convert from VCF format to tab-delimited text file <br>
zcat file.vcf.gz | <a href="perl_module.html#vcf-to-tab">vcf-to-tab</a> &gt; out.tab <br>
</p>
</div>


<div class="item">
<a name="filter" class="Q">Filtering</a>
<p></p>
<p class="codebox">
# Filter by QUAL and minimum depth <br>
<a href="perl_module.html#vcf-annotate">vcf-annotate</a> --filter Qual=10/MinDP=20 
</p>
</div>


<div class="item">
<a name="isec" class="Q">Intersections, complements</a>
<p></p>
<p class="codebox">
# Include positions which appear in at least two files  <br>
<a href="perl_module.html#vcf-isec">vcf-isec</a> -o -n +2 A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c &gt; out.vcf.gz <br>
<br>
# Exclude from A positions which appear in B and/or C   <br>
<a href="perl_module.html#vcf-isec">vcf-isec</a> -c A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c &gt; out.vcf.gz   <br>
<br>
# Fast hstlib implementation
<a href="htslib.html#isec">vcf isec</a> -n =2 A.vcf.gz B.vcf.gz <br>
</p>
</div>


<div class="item">
<a name="merge" class="Q">Merging</a>
<p></p>
<p class="codebox">
<a href="perl_module.html#vcf-merge">vcf-merge</a> A.vcf.gz B.vcf.gz | bgzip -c &gt; C.vcf.gz <br>
<a href="htslib.html#merge">vcf merge</a> A.vcf.gz B.vcf.gz <br>
</p>
</div>

<div class="item">
<a name="isec" class="Q">Querying</a>
<p></p>
<p class="codebox"><a href="perl_module.html#vcf-query">vcf-query</a> file.vcf.gz 1:10327-10330 -c NA0001
</p>
</div>



<div class="item">
<a name="shuffle" class="Q">Reordering columns</a>
<p></p>
<p class="codebox"><a href="perl_module.html#vcf-shuffle-cols">vcf-shuffle-cols</a> -t template.vcf.gz file.vcf.gz > out.vcf
</p>
</div>


<div class="item">
<a name="stats" class="Q">Stats</a>
<p></p>
<p class="codebox"><a href="perl_module.html#vcf-stats">vcf-stats</a> file.vcf.gz <br>
<a href="htslib.html#check">vcf check</a> file.vcf.gz &gt; file.vchk && plot-vcfcheck file.vchk -p plot/<br>
</p>
</div>


<div class="item">
<a name="subset" class="Q">Stripping columns</a>
<p></p>
<p class="codebox"><a href="perl_module.html#vcf-subset">vcf-subset</a> -c NA0001,NA0002 file.vcf.gz | bgzip -c &gt; out.vcf.gz
</p>
</div>



<div class="item">
<a name="one-liners" class="Q">Useful shell one-liners</a>
<p> This sections lists some usefull one line commands. Note that there are also dedicated convenience scripts
<a href="perl_module.html#vcf-sort">vcf-sort</a> and <a href="perl_module.html#vcf-concat">vcf-concat</a> 
which do the same but also perform some basic sanity checks.
All examples in BASH.</p>
<p class="codebox">
# Replace VCF header. The file must be compressed by bgzip. <br>
tabix -r header.txt in.vcf.gz > out.vcf.gz                  <br>
<br>
# Sort VCF file keeping the header. The head command is for performance.    <br>
(zcat file.vcf.gz | head -100 | grep ^#;                                    <br>
 zcat file.vcf.gz | grep -v ^# | sort -k1,1d -k2,2n;)  \                    <br>
 | bgzip -c &gt; out.vcf.gz <br>
<br>
# Merge (that is, concatenate) two VCF files into one, keeping the header   <br>
#   from first one only.                    <br>
(zcat A.vcf.gz | head -100 | grep ^#; \     <br>
 zcat A.vcf.gz | grep -v ^#;          \     <br>
 zcat B.vcf.gz | grep -v ^#; )        \     <br>
 | bgzip -c &gt; out.vcf.gz     
</p>
</div>


<div class="item">
<a name="validator" class="Q">VCF validation</a>
<p> Both <span class="cmd">vcftools</span> and <span class="cmd">Vcf.pm</span> can be used for validation. The
    first validates VCFv4.0, the latter is able to validate the older versions as well.
</p>
<p class="codebox">
perl -M<a href="perl_module.html#Vcf.pm">Vcf</a> -e validate example.vcf  <br>
perl -I/path/to/the/module/ -M<a href="perl_module.html#Vcf.pm">Vcf</a> -e validate example.vcf   <br>
<a href="perl_module.html#vcf-validator">vcf-validator</a> example.vcf
</p>
</div>


<div class="item">
<a name="more" class="Q">...and more</a>
<p>This page gives just a list of basic capabilities. For more, please go to the <a href="options.html">vcftools's options</a> page
and the <a href="perl_module.html">Perl API and scripts</a> page.</p>
</div>

