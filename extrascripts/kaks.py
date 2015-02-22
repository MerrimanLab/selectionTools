import vcf
import sys

# HOW TO RUN:
# use snpEff to annoate VCF file
# SLOW:
# python kaks.py annotated_vcf.vcf > kaks.txt
# FASTER:
# grep '^#\|missense_variant\|synonymous_variant' annotated_vcf.vcf > mis_syn.txt 
# python kaks.py mis_syn.txt > kaks.txt

file = sys.argv[1]

vcf_reader = vcf.Reader(open(file, 'r'))
gene_name_index = vcf_reader.infos['ANN'][3].split("|").index(" Gene_Name ")
geneid_index =vcf_reader.infos['ANN'][3].split("|").index(" Gene_ID ")
annotation_index = vcf_reader.infos['ANN'][3].split("|").index(" Annotation ")
feature_id_index = vcf_reader.infos['ANN'][3].split("|").index(" Feature_ID ")

#key = geneid, value = [ka, ks, gene_name]
genes_ka_ks = {}

def add_gene(record):
    curr_gene = str(record.INFO['ANN']).split('|')[geneid_index]
    for r in record.INFO['ANN']:
        gene = str(r).split('|')[geneid_index]
        #print gene
        #print str(r).split('|')[annotation_index]
        if gene in genes_ka_ks:
            if str(r).split('|')[annotation_index] == "synonymous_variant" :
                genes_ka_ks[gene] = [genes_ka_ks[gene][0], genes_ka_ks[gene][1] + record.num_het + 2 * record.num_hom_alt, 
                                     genes_ka_ks[gene][2]]
            
            if str(r).split('|')[annotation_index] == "missense_variant" :
                genes_ka_ks[gene] = [genes_ka_ks[gene][0] + record.num_het + 2 * record.num_hom_alt, genes_ka_ks[gene][1],
                                     genes_ka_ks[gene][2]]
        
        else:
            if str(r).split('|')[annotation_index] == "synonymous_variant" :
                genes_ka_ks[gene] = [0,  record.num_het + 2 * record.num_hom_alt, 
                                     str(r).split('|')[gene_name_index]]
            
            if str(r).split('|')[annotation_index] == "missense_variant" :
                genes_ka_ks[gene] = [record.num_het + 2 * record.num_hom_alt, 0, 
                                     str(r).split('|')[gene_name_index]]


print "GeneID\tGeneName\tka\tks\tka_div_ks_plus1"
for record in vcf_reader:
    add_gene(record)
    
for gene in genes_ka_ks:
    print( gene +"\t"+ genes_ka_ks[gene][2] +"\t"+ str(genes_ka_ks[gene][0]) 
          +"\t"+ str(genes_ka_ks[gene][1]) +"\t"+ str(genes_ka_ks[gene][0]/float(genes_ka_ks[gene][1] +1)) ) 


            
        
