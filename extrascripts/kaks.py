import vcf
import sys

# HOW TO RUN:
# use snpEff to annoate VCF file
# grep '^#\|missense_variant\|synonymous_variant' annotated_vcf.vcf > mis_syn.txt 
# python kaks.py mis_syn.txt > kaks.txt



file = sys.argv[1]
vcf_reader = vcf.Reader(open(file, 'r'))
gene_name_index = vcf_reader.infos['ANN'][3].split("|").index(" Gene_Name ")
geneid_index =vcf_reader.infos['ANN'][3].split("|").index(" Gene_ID ")
annotation_index = vcf_reader.infos['ANN'][3].split("|").index(" Annotation ")
feature_id_index = vcf_reader.infos['ANN'][3].split("|").index(" Feature_ID ")

#ka = number alleles with missense variants
ka = 0
# ks = number of alleles with synonymous variants
ks = 0
# total number of alleles for variant
total = 0


#record = vcf_reader.next()
#while (int(record.POS) < 24852068) :
i = 1
print "GeneID\tGeneName\tka\tks\tka_div_ks_plus1"
for record in vcf_reader:
    curr_gene = str(record.INFO['ANN']).split('|')[geneid_index]
    if i == 1:
        prev_gene = curr_gene
        i += 1    
    if curr_gene != prev_gene:
        print curr_gene + "\t" + str(record.INFO['ANN']).split('|')[gene_name_index] +"\t" + str(ka) +"\t" + str(ks) +"\t" + str( (ka) / float(ks+1))
        #print "total " + str(total)
        #make sure to reset ka, ks and total at each gene change
        ka =0; ks =1; total =0
    prev_gene = curr_gene
    for r in record.INFO['ANN']:
        #print( str(r).split('|')[geneid_index] +" " + str(r).split('|')[feature_id_index] +
        #      " " + str(r).split('|')[annotation_index] )
        if str(r).split('|')[annotation_index] == "synonymous_variant" :
            #print(str(record.num_called) + " " + str(record.num_hom_ref) + " " + str(record.num_het) + " " +str(record.num_hom_alt))
            ks += record.num_het + 2 * record.num_hom_alt
            total += 2 * record.num_called
            #print ks
        if str(r).split('|')[annotation_index] == "missense_variant" :
            #print(str(record.num_called) + " " + str(record.num_hom_ref) + " " + str(record.num_het) + " " +str(record.num_hom_alt))
            ka += record.num_het + 2 * record.num_hom_alt
            total += 2 * record.num_called
            #print ks
            
        
