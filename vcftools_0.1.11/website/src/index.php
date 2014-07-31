<?php

function validate_path($path)
{
    if ( preg_match("/^[A-Za-z0-9_]+$/",$path) && file_exists("src/$path.inc") ) return $path;
    return "index";
}

$titles = array(
        'index'         => 'VCFtools',
        'perl_module'   => 'VCFtools: Perl tools and API',
        'htslib'        => 'VCFtools: htslib VCF commands',
        'docs'          => 'VCFtools Documentation',
        'license'       => 'VCFtools License',
        'specs'         => 'VCF Specification',
        'links'         => 'VCF Links',
        'options'       => 'vcftools Options',
        );

if (isset($argc)) { $_GET['pg']=$argv[1]; }
$path  = array_key_exists('pg',$_GET) ? validate_path($_GET['pg']) : 'index';
$title = array_key_exists($path,$titles) ? $titles[$path] : $titles['index'];


?>
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
    <script type="text/javascript">
        var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
        document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
    </script>
    <script type="text/javascript">
        try {
            var pageTracker = _gat._getTracker("UA-272183-4");
            pageTracker._trackPageview();
        } catch(err) {}
    </script>

    <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
    <link rel="stylesheet" type="text/css" href="default.css" media="screen">
    <link href='favicon.png' rel='shortcut icon' type='image/png'>
    <link href='favicon.png' rel='icon' type='image/png'>
    
    <title><?php echo $title; ?></title>
</head>



<body>
<div class="container">
<div class="main">
    <div class="header">
    <div class="title">
        <a href="index.html">VCFtools</a>
    </div>
    </div>
<div class="content">

<?php 
    include("$path.inc");
?>

</div>

<div class="sidenav">
    <h1>Navigation</h1>
    <ul>
    <li><a href="index.html">Main</a></li>
    <li><a href="http://sourceforge.net/projects/vcftools/">Sourceforge page</a></li>
    <li><a href="docs.html">Documentation</a></li>
    <li><a href="license.html">License</a></li>
    <li><a href="specs.html">VCF specification</a></li>
    <li><a href="links.html">Links</a></li>
    <li><a href="http://www.1000genomes.org/">1000 Genomes</a></li>
    </ul>
</div>
<div class="clearer"><span></span></div>
</div>
</div>
</body></html>
