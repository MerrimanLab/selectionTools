
<div class="item">
<h1>The htslib VCF commands</h1>
<p> <a href="https://github.com/samtools/htslib">HTSlib</a> is a C library for
high-throughput sequencing data formats. It is designed for speed and works with both VCF and 
<A href="http://www.1000genomes.org/wiki/analysis/variant-call-format/bcf-binary-vcf-version-2">BCFv2</A>.
</p>

<h3>Download and installation</h3>
The library is hosted on github. It can be downloaded and compiled the usual way:
<p class="codebox">
     # Download. The first command is run only once. The second
<br> # is run whenever the latest snapshot from github is needed.
<br> git clone git://github.com/samtools/htslib.git htslib
<br> git pull
<br>
<br> # Compile
<br> cd htslib
<br> make git-stamp
<br>
<br> # Run
<br> htscmd vcfcheck file.vcf.gz
<br>
<br> # Symlinking the executable as 'vcf' allows terse execution of VCF commands:
<br> ln -s htscmd vcf
<br> vcf check file.vcf.gz
</p>

<h3>The tools</h3>
<ul class="faq">
<li> <a href="#check">vcf check</a>
<li> <a href="#filter">vcf filter</a>
<li> <a href="#gtcheck">vcf gtcheck</a>
<li> <a href="#isec">vcf isec</a>
<li> <a href="#merge">vcf merge</a>
<li> <a href="#norm">vcf norm</a>
<li> <a href="#query">vcf query</a>
<li> <a href="#view">vcf view</a>
</ul>

</div>

<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js"></script>
<script type="text/javascript">
jQuery(document).ready(function() {
        jQuery(".usageText").hide();
        jQuery(".usageToggle").click(function() { jQuery(this).next(".usageText").slideToggle(100); });
    });
</script>

<div class="item">
<a name="check" class="Q">vcf check</a>
<p> Extract stats from a VCF/BCF file or compare two VCF/BCF files. The resulting text file can be plotted using 
<span class="cmd">plot-vcfcheck</span>
</p>
<p class="codebox">
vcf check file.vcf.gz &gt; file.vchk 
<br> plot-vcfcheck file.vchk -p plots/</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("vcf check 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="filter" class="Q">vcf filter</a>
<p> Filtering by both supervised and unsupervised learning.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("vcf filter 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="gtcheck" class="Q">vcf gtcheck</a>
<p> A tool for detecting sample swaps and contaminations
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("vcf gtcheck 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="isec" class="Q">vcf isec</a>
<p> Fast alternative to <A href="perl_module.html#vcf-isec">vcf-isec</A>
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("vcf isec 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="merge" class="Q">vcf merge</a>
<p> Fast alternative to <A href="perl_module.html#vcf-merge">vcf-merge</A>
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("vcf merge 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="norm" class="Q">vcf norm</a>
<p> Left-align and normalize indels to the shortest possible representation.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("vcf norm 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="query" class="Q">vcf query</a>
<p> Fast alternative to <A href="perl_module.html#vcf-query">vcf-query</A>
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("vcf query 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>


