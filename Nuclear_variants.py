#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import argparse
import yaml
import collections
import itertools
import xlwt

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqIO import parse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():

    # set groups and samples dict
    # samples[sample] = group
    # groups[group] = sample_list
    groups = dict()
    samples = dict()
    for group in config['pops']:
        groups[group] = config['pops'][group]
        for sample in config['pops'][group]:
            samples[sample] = group
    group_list = sorted(groups)
    skip_samples = set(config['skip_samples'])
    
    # set all comparisons between 2 pops
    cmps = [cmp for cmp in itertools.combinations(group_list, 2)]
    print cmps

    if args.verbose :
        print "skip_samples", skip_samples
        print "samples", samples
        print "groups", groups
        print "cmps", cmps

    d1 = "{}/{}".format(config['top_dir'], config['res_dir'])
    if not os.path.isdir(d1):
        os.makedirs(d1)

    vcf_1 = "{}/{}/F1.vcf".format(config['top_dir'], config['res_dir'])

    data = dict()

    # read vcf and count alleles per group
    # data[contig]['sites'][pos]['ref_alleles'] : list
    # data[contig]['sites'][pos]['alt_alleles'] : list
    # data[contig]['sites'][pos][group + '_alleles'] : alleles_set (nucl)
    # data[contig]['sites'][pos][group][sample] : genotype string (ex "A/C")
    read_vcf(data, vcf_1, groups, samples, skip_samples)
    exit(0)
    #print "{:,} contigs".format(len(data))

    # annot variants with groups comparison
    # data[contig]['sites'][pos][cmp1] = O/1          1 = site is selected for this cmp
    # data[contig]['sites'][pos]['cmp_select'] = n    n = sum of selected cmp per site (with at least one selected cmp)
    # data[contig]['sites_select'] = n                n = sum of selected sites per contig (with at least one selected site)
    # data[contig][cmp1]['sites_select'] = n          n = sum of selected sites per contig and cmp
    annot_variants_cmp(data, cmps)

    # annot contigs with TAIR
    # data[contig]['TAIR'] = dict, keys=id, symbol, description, genome_pos, strand, length
    annot_contigs(data)

    # count results
    count = collections.defaultdict(int)
    selected_TAIR = set()
    selected_contigs = set()
    for contig in data:
        count['contigs'] +=1
        #print "{} {} / {} sites".format(contig,  data[contig]['select'], len(data[contig]['sites']))
        count['sites'] += len(data[contig]['sites'])
        if data[contig]['sites_select'] > 0:
            count['s_sites'] += data[contig]['sites_select']
            count['s_contigs'] += 1
            selected_contigs.add(contig)
            if 'TAIR' in data[contig]:
                count['TAIR'] +=1
                selected_TAIR.add(data[contig]['TAIR']['id'])

    print "{:,} polymorphic sites in {:,} contigs".format(count['sites'], count['contigs'])
    print "select : {:,} sites in {:,} contigs {:,} contigs with TAIR annot {:,} TAIR genes".format(count['s_sites'], count['s_contigs'], count['TAIR'], len(selected_TAIR))

    # annot selected variants with codon syn / non-syn (search codon in ref fasta file)
    # data[contig]['fasta'] : fasta record
    # data[contig]['sites'][pos]['ref_codons'] : list
    # data[contig]['sites'][pos]['alt_codons'] : list
    # data[contig]['sites'][pos]['ref_aa']     : list
    # data[contig]['sites'][pos]['alt_aa']     : list
    # data[contig]['sites'][pos]['syn']        : list
    annot_codons(data)
    
    # create excel file
    print "write", args.xls
    wb = xlwt.Workbook()
    ws1 = wb.add_sheet("contigs details")
    ws2 = wb.add_sheet("polymorphic sites details")
    fmt_title = xlwt.easyxf('font: name Times New Roman, color-index black, bold on')

    TAIR_keys = ['id', 'symbol', 'description', 'genome_pos', 'strand', 'length']

    ############## sheet1 (per contig) ##############
    # titles sheet 1 (per contig)
    l1 = 0; c1 = 0
    cols = ["contig",]
    for k in TAIR_keys:
        cols.append("TAIR {}".format(k))
    cols.append("# polymorphic sites")
    cols.append("# selected sites")
    for cmp1 in cmps:
        G1, G2 = cmp1
        cols.append("{}-{} ssites".format(G1, G2)) # selected sites per cmp
    for col in cols:
        ws1.write(l1, c1, col, fmt_title); c1 += 1
    l1 += 1; c1 = 0

    # lines sheet 1 (per contig)
    for contig in sorted(data):
        # skip unselected contigs
        if data[contig]['sites_select'] == 0:
            continue

        # collect values
        cols = [contig, ]
        if 'TAIR' in data[contig]:
            for k in TAIR_keys:
                cols.append(str(data[contig]['TAIR'][k]))
        else:
            for k in TAIR_keys:
                cols.append("-")
        cols.append(len(data[contig]['sites'])) # polymorphic sites
        cols.append(data[contig]['sites_select']) # total selected sites 
        for cmp1 in cmps:
            cols.append(data[contig][cmp1]['sites_select']) # selected sites per cmp

        # print line
        for col in cols:
            ws1.write(l1, c1, col); c1 += 1
        l1 += 1; c1 = 0

    ############## sheet2 (per site) ##############
    group_list = sorted(list(groups))
    # titles sheet 2 (per site)
    l2 = 0; c2 = 0
    cols = ["contig", "pos"]
    for k in TAIR_keys:
        cols.append("TAIR {}".format(k))
    cols.append("Ref alleles")
    cols.append("Ref codons")
    cols.append("Ref AAs")
    cols.append("Alt alleles")
    cols.append("Atl codons")
    cols.append("Atl AAs")
    cols.append("S/NS")
    for group in group_list:
        cols.append("{} allele".format(group))
    cols.append("# cmp selected") # nb of select cmp per site
    for cmp1 in cmps:
        G1, G2 = cmp1
        cols.append("{}-{}".format(G1, G2)) # cmp1 select 0/1
    for col in cols:
        ws2.write(l2, c2, col, fmt_title); c2 += 1
    l2 += 1; c2 = 0

    # lines sheet 2 (per site)
    for contig in sorted(data):
        # skip unselected contigs
        if data[contig]['sites_select'] == 0:
            continue
        
        for pos in sorted(list(data[contig]['sites'])):
            # skip unselected sites
            if data[contig]['sites'][pos]['cmp_select'] == 0:
                continue
            
            # collect values
            cols = [contig, pos]
            if 'TAIR' in data[contig]:
                for k in TAIR_keys:
                    cols.append(str(data[contig]['TAIR'][k]))
            else:
                for k in TAIR_keys:
                    cols.append("-")
            keys = ['ref_alleles', 'ref_codons', 'ref_aa', 'alt_alleles', 'alt_codons', 'alt_aa', 'syn']
            for key in keys:
                cols.append(",".join(data[contig]['sites'][pos][key]))
            for group in group_list:
                cols.append(",".join(sorted(list(data[contig]['sites'][pos][group + '_alleles']))))
            cols.append(data[contig]['sites'][pos]['cmp_select']) # nb of select cmp per site
            for cmp1 in cmps:
                cols.append(data[contig]['sites'][pos][cmp1]) # cmp1 select 0/1

            # print line
            for col in cols:
                ws2.write(l2, c2, col); c2 += 1
            l2 += 1; c2 = 0

    # write excel file
    wb.save(args.xls)

def read_vcf(data, vcf, groups, samples, skip_samples):
    # read vcf and count alleles per group
    # data[contig]['sites'][pos]['ref_alleles'] : list
    # data[contig]['sites'][pos]['alt_alleles'] : list
    # data[contig]['sites'][pos][group + '_alleles'] : alleles_set (nucl)
    # data[contig]['sites'][pos][group][sample] : genotype string (ex "A/C")

    sample2col = dict()
    group_list = sorted(list(groups))
    
    count_ref_alleles_size = collections.Counter()
    count_alt_alleles_size = collections.Counter()
    count_ref_alleles_nt = collections.Counter()
    count_alt_alleles_nt = collections.Counter()
    with open(vcf) as f_in:
        for l in f_in:
            if l.startswith('#CHROM'):
                tab = l.strip().split('\t')
                for i in range(9, len(tab)):
                    sample = tab[i].replace('.bam', '')
                    if sample in skip_samples:
                        continue
                    #group = samples[sample]
                    sample2col[sample] = i
                    #print "{} {} {}".format(i, sample, group)
                    #group_cols[group].append(i)
                #print "group_cols", group_cols
                continue

            if l.startswith('#'):
                continue

            tab = l.strip().split('\t')

            contig = tab[0]
            pos = int(tab[1])
            #if contig == 'c26585_g1_i1|m.23326' and pos == 91:
            #    print "GOT IT", pos
            #else:
            #    continue
            if contig not in data:
                data[contig] = dict()
                data[contig]['sites'] = dict()
                #if len(data) >= args.max_contigs:
                #    break
            if pos not in data[contig]['sites']:
                data[contig]['sites'][pos] = dict()
                # collect alleles for this variant
                ref_alleles = tab[3].split(',')
                alt_alleles = tab[4].split(',')
                count_ref_alleles_size[len(ref_alleles)] += 1
                count_alt_alleles_size[len(alt_alleles)] += 1
                for a in ref_alleles:
                    count_ref_alleles_nt[a] += 1
                for a in alt_alleles:
                    count_alt_alleles_nt[a] += 1
                data[contig]['sites'][pos]['ref_alleles'] = ref_alleles
                data[contig]['sites'][pos]['alt_alleles'] = alt_alleles

            # collect alleles per group
            for group in group_list:
                data[contig]['sites'][pos][group] = dict()  # per sample genotypes
                data[contig]['sites'][pos][group + "_alleles"] = set() # all alleles for this group
                for sample in groups[group]:
                    i = sample2col[sample]
                    geno = tab[i]
                    if geno != ".|.":
                        A1, A2 = alleles_recode(contig, pos, geno, ref_alleles, alt_alleles, l)
                        data[contig]['sites'][pos][group][sample] = "{}/{}".format(A1, A2)
                        data[contig]['sites'][pos][group + "_alleles"] |= set([A1, A2])
            #print "contig {} pos {}".format(contig, pos)
            #for group in group_list:
            #    print "\tgroup {}, alleles {}, {} genotypes {}".format(group, data[contig]['sites'][pos][group + "_alleles"], len(data[contig]['sites'][pos][group]), data[contig]['sites'][pos][group])
            #exit(0)
            #print data[contig]


    print "count_ref_alleles_size", count_ref_alleles_size
    print "count_alt_alleles_size", count_alt_alleles_size
    print "count_ref_alleles_nt", count_ref_alleles_nt
    print "count_alt_alleles_nt", count_alt_alleles_nt
    return data

def alleles_recode(contig, pos, geno, ref_allele, alt_alleles, l):
    # !!! big hypothesis : genotypes are numbererd in the same order as alt alleles, 
    # example ref_allele A, alt_alleles C,G : geno 0|1 = A/C, geno 0|2 = A/G
    # !!! ref alleles may be '-'
    recode = { '0' : ref_allele[0] }
    for i in range(len(alt_alleles)):
        recode[str(i + 1)] = alt_alleles[i]

    a1, a2 = geno.split('|')
    A1 = recode[a1]
    A2 = recode[a2]
    #if len(alt_alleles) > 3:
    #    print "{} {} ref_allele {}, alt_alleles {}, recode {}, {}/{} => {}/{}".format(contig, pos, ref_allele, alt_alleles, recode, a1, a2, A1, A2)
    return A1, A2

def annot_variants_cmp(data, cmps):
    # data[contig]['sites'][pos][cmp1] = O/1          1 = site is selected for this cmp
    # data[contig]['sites'][pos]['cmp_select'] = n    n = sum of selected cmp per site (with at least one selected cmp)
    # data[contig]['sites_select'] = n                n = sum of selected sites per contig (with at least one selected site)
    # data[contig][cmp1]['sites_select'] = n          n = sum of selected sites per contig and cmp

    for contig in data:
        data[contig]['sites_select'] = 0
        for cmp1 in cmps:
            data[contig][cmp1] = {'sites_select' : 0}
        for pos in data[contig]['sites']:
            data[contig]['sites'][pos]['cmp_select'] = 0
            for cmp1 in cmps:
                data[contig]['sites'][pos][cmp1] = 0
                G1, G2 = cmp1
                
                # check min_geno_per_group
                if len(data[contig]['sites'][pos][G1]) < config['min_geno_per_group']:
                    continue
                if len(data[contig]['sites'][pos][G2]) < config['min_geno_per_group']:
                    continue

                # check fixed in G1                
                if len(data[contig]['sites'][pos][G1 + "_alleles"]) != 1:
                    continue
                A1 = list(data[contig]['sites'][pos][G1 + "_alleles"])[0]

                # check fixed in G2
                if len(data[contig]['sites'][pos][G2 + "_alleles"]) != 1:
                    continue
                A2 = list(data[contig]['sites'][pos][G2 + "_alleles"])[0]

                # check diff between groups
                if A1 == A2:
                    continue
                
                # mark selected sites for this cmp
                data[contig]['sites'][pos][cmp1] = 1
                # count selected cmp per site
                data[contig]['sites'][pos]['cmp_select'] += 1
                # count selected sites per cmp
                data[contig][cmp1]['sites_select'] += 1
                #print "contig {} pos {} cmp {} selected A1={} A2={}".format(contig, pos, cmp1, A1, A2)

            # count selected sites per contig (with at least one selected site)
            if data[contig]['sites'][pos]['cmp_select'] > 0:
                data[contig]['sites_select'] += 1
                # print "contig {} pos {} nb of selection {}".format(contig, pos, data[contig]['sites'][pos]['cmp_select'])

def annot_contigs(data):
    """ 
    TAIR_annots example
    c11248_g1_i1|m.7237	AT3G53710.2	57.73	485	143	10	10	1362	3	459	4e-167	 483
    CONTIG              TAIR ID

    TAIR_details example
    >AT1G51370.2 | Symbols:  | F-box/RNI-like/FBD-like domains-containing protein | chr1:19045615-19046748 FORWARD LENGTH=346
    TAIR_ID        SYMBOL      DESCRIPTION                                          GENOME_POS             STRAND  LENGTH

    # data[contig]['TAIR'] = dict, keys=id, symbol, description, genome_pos, strand, length
    """

    TAIR_annots = config['TAIR_annots']
    TAIR_details = config['TAIR_details']

    # /gepv/projects1/users/lenou/LN-data/Analyses_Nutans/refE/GO_enrichment/README
    # nohup blastx -db db/At_pep -query ../ref/Assem.ORF.NutansE.fasta -max_target_seqs 1 -evalue 1e-5 -outfmt 6 -out blast_out/blast_ref_TAIR.out_sfi &
    # blastx : query=contig, target=tair, max target per query = 1 
    # => per contig : 1 tair_id
    # => per tair_id : possibly n contigs

    contig2tair_id = dict()
    tair_id2contigs = collections.defaultdict(list)
    #print "read", TAIR_annots
    with open(TAIR_annots) as f_in:
        for l in f_in:
            tab = l.strip().split('\t')
            contig = tab[0]
            tair_id = tab[1]
            if contig in data:
                #print "tair_id {} => contig {}".format(tair_id, contig)
                data[contig]['TAIR'] = dict()
                data[contig]['TAIR']['id'] = tair_id
                tair_id2contigs[tair_id].append(contig)
                contig2tair_id[contig] = tair_id

    #print "{} contigs with a tair_id".format(len(contig2tair_id))
    #print "{} tair_id with a contig".format(len(tair_id2contigs))
    #for tair_id in tair_id2contigs:
    #    print "tair_id {} => contigs {}".format(tair_id, tair_id2contigs[tair_id])
    #    if len(tair_id2contigs[tair_id]) > 1:
    #        print "tair_id {} > 1 contigs {}".format(tair_id, tair_id2contigs[tair_id])

    #print "read", TAIR_details
    with open(TAIR_details) as f_in:
        for l in f_in:
            tair_id, symbol, description, tmp = l.strip().split('|')
            # strip() to remove first and last space chars
            # tair_id : first char '>'
            tair_id = tair_id.replace('>', '').strip()
            symbol = symbol.replace('Symbols:', '').strip()
            description = description.strip()
            # there is a space tab at the begenning ==> add a dummy var
            dummy, genome_pos, strand, length = tmp.split(' ')
            length = int(length.replace('LENGTH=', ''))
            if tair_id in tair_id2contigs:
                #print "!!1 tair_id {} => contig {}".format(tair_id, tair_id2contigs[tair_id])
                for contig in tair_id2contigs[tair_id]:
                    data[contig]['TAIR']['symbol'] = symbol
                    data[contig]['TAIR']['description'] = description
                    data[contig]['TAIR']['genome_pos'] = genome_pos
                    data[contig]['TAIR']['strand'] = strand
                    data[contig]['TAIR']['length'] = length
                    #print data[contig]['TAIR']
            
    # check
    count = collections.defaultdict(int)
    for contig in data:
        if 'TAIR' in data[contig]:
            count['TAIR_annots'] +=1
            if 'description' in data[contig]['TAIR']:
                count['TAIR_details'] +=1
    #print "{:,} contigs, {:,} with TAIR annot, {:,} with TAIR details".format(len(data), count['TAIR_annots'], count['TAIR_details'])


def annot_codons(data):
    """ annot selected variants with codon syn / non-syn (search codon in ref fasta file) """
    # data[contig]['fasta'] : fasta record
    # data[contig]['sites'][pos]['ref_codons'] : list
    # data[contig]['sites'][pos]['alt_codons'] : list
    # data[contig]['sites'][pos]['ref_aa']     : list
    # data[contig]['sites'][pos]['alt_aa']     : list
    # data[contig]['sites'][pos]['syn']        : list
    
    # set of selected contigs
    selected_contigs = set()
    for contig in sorted(data):
        if data[contig]['sites_select'] > 0:
            selected_contigs.add(contig)
    print "{} selected contigs".format(len(selected_contigs))

    # collect fasta for selected contigs
    in_fasta = config['ref']
    n_in = 0
    n_out = 0
    print "read", in_fasta
    with open(in_fasta) as f_in:
        for record in SeqIO.parse(f_in, "fasta") :
            n_in += 1
            if record.id in selected_contigs:
                data[record.id]['fasta'] = record
                n_out += 1
            fasta_file = "{}/{}/{}.fasta".format(config['top_dir'], config['res_dir'], record.id.replace('|', '_'))
            print "write", fasta_file
            SeqIO.write(record, fasta_file, "fasta")
    print "{}/{} contigs from {}".format(n_out, n_in, in_fasta)

    # collect annotations for selected sites
    for contig in sorted(selected_contigs):
        if 'fasta' not in data[contig]:
            print "!! missing fasta for contig", contig
            exit(1)
        ref_seq = data[contig]['fasta']
        for pos in sorted(list(data[contig]['sites'])):
            # skip unselected sites
            if data[contig]['sites'][pos]['cmp_select'] == 0:
                continue
            #print "work on", contig, pos
            vcf_ref_alleles = data[contig]['sites'][pos]['ref_alleles']
            vcf_alt_alleles = data[contig]['sites'][pos]['alt_alleles']
            # cf /gepv/pgen/silene_cg_cyto/analyse_2019-01/5c_variants_NS_v1.py
            pos0 = pos - 1             # !! reads2snp vcf pos : 1-based, BioPython seq pos : 0-based
            # check
            fasta_ref_allele = ref_seq[pos0]
            if vcf_ref_alleles[0] != fasta_ref_allele:
                print "!! ERROR contig {} pos {} pos0 {} phase {} vcf ref_allele {} fasta ref_allele {}".format(contig, pos, pos0, phase, vcf_ref_alleles[0], fasta_ref_allele)
                exit(1)
            
            phase = pos0 % 3

            data[contig]['sites'][pos]['ref_codons'] = list()
            data[contig]['sites'][pos]['ref_aa'] = list()
            for vcf_ref_allele in vcf_ref_alleles:
                if phase == 0:
                    codon = vcf_ref_allele    + ref_seq[pos0 + 1] + ref_seq[pos0 + 2]
                elif phase == 1:
                    codon = ref_seq[pos0 - 1] + vcf_ref_allele    + ref_seq[pos0 + 1]
                elif phase == 2:
                    codon = ref_seq[pos0 - 2] + ref_seq[pos0 - 1] + vcf_ref_allele
                if '-' in codon:
                    aa = "-"
                else:
                    aa = str(Seq(codon).translate(table=config['NCBI_trans_table']))
                data[contig]['sites'][pos]['ref_codons'].append(codon)
                data[contig]['sites'][pos]['ref_aa'].append(aa)

            data[contig]['sites'][pos]['alt_codons'] = list()
            data[contig]['sites'][pos]['alt_aa'] = list()
            data[contig]['sites'][pos]['syn'] = list()
            for vcf_alt_allele in vcf_alt_alleles:
                if phase == 0:
                    codon  = vcf_alt_allele   + ref_seq[pos0 + 1] + ref_seq[pos0 + 2]
                elif phase == 1:
                    codon = ref_seq[pos0 - 1] + vcf_alt_allele    + ref_seq[pos0 + 1]
                elif phase == 2:
                    codon = ref_seq[pos0 - 2] + ref_seq[pos0 - 1] + vcf_alt_allele
                if '-' in codon:
                    aa = "-"
                else:
                    aa = str(Seq(codon).translate(table=config['NCBI_trans_table']))
                data[contig]['sites'][pos]['alt_codons'].append(codon)
                data[contig]['sites'][pos]['alt_aa'].append(aa)
                if data[contig]['sites'][pos]['alt_aa'][0] == data[contig]['sites'][pos]['ref_aa'][0]:
                    data[contig]['sites'][pos]['syn'].append('S')
                else:
                    data[contig]['sites'][pos]['syn'].append('NS')

            keys = ['ref_alleles', 'ref_codons', 'ref_aa', 'alt_alleles', 'alt_codons', 'alt_aa', 'syn']
            tmp = [ "{}={}".format(i, data[contig]['sites'][pos][i]) for i in keys]
            #print "contig {} pos {} pos0 {} phase {} {}".format(contig, pos, pos0, phase, ", ".join(tmp))

if __name__ == '__main__':
    description = """ 
    Select variants from reads2snp vcf result
    for each comparison between 2 pops, select variant fixed in pop and diff between pop
    Determine synonymous or non synonymous variants type
    annot with tair gene and description
    Select variants seggregation for each comparison between 2 pops
    Create xls file with 2 sheets for contigs details and polymorphic sites details
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--verbose", action="store_const", const="-v", default="")
    parser.add_argument("--config", help="yaml config file")
    parser.add_argument("--xls", help="output xls file", default="nucl_variants.xls")
    parser.add_argument("--max_contigs", type=int, default="10")
    
    args = parser.parse_args()
    config = yaml.load(open(args.config))

    main() 
