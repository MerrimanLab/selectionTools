import vcf
file = "mis_syn.txt"
vcf_reader = vcf.Reader(open(file, 'r'))
gene_name_index = vcf_reader.infos['ANN'][3].split("|").index(" Gene_Name ")
geneid_index =vcf_reader.infos['ANN'][3].split("|").index(" Gene_ID ")
annotation_index = vcf_reader.infos['ANN'][3].split("|").index(" Annotation ")
feature_id_index = vcf_reader.infos['ANN'][3].split("|").index(" Feature_ID ")

ka = 0
ks = 0
total = 0
record = vcf_reader.next()
while (int(record.POS) < 24852068) :
    curr_gene = str(record.INFO['ANN']).split('|')[geneid_index]
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
            
    record= vcf_reader.next()
    next_gene = str(record.INFO['ANN']).split('|')[geneid_index]
    if curr_gene != next_gene:
        print curr_gene + " " + str(record.INFO['ANN']).split('|')[gene_name_index]
        print "ks " + str(ks)
        print "ka " + str(ka)
        print "total " + str(total)
        ka =0; ks =0; total =0
