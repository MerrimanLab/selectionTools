
FILES=docs index license options perl_module specs links htslib

install:
	    @for i in $(FILES); do echo "php src/index.php $$i > $$i.html"; php src/index.php $$i > $$i.html; done;

