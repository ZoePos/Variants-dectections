#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import collections
import itertools
import xlwt

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqIO import parse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():

    # read sample and set group_list (ie population list)
    samples, groups = read_samples(args.samples)
    group_list = sorted(groups.keys())
    group_list.remove('all')
    group_list.remove('W')

    # set all comparisons between 2 pops
    cmps = [cmp for cmp in itertools.combinations(group_list, 2)]


    # create excel file with 2 sheets and formats for titles
    wb = xlwt.Workbook()
    ws1 = wb.add_sheet("summary count")
    ws2 = wb.add_sheet("polymorphic sites details")
    fmt_title = xlwt.easyxf('font: name Times New Roman, color-index black, bold on')

    # write titles to sheet1
    l1 = 0; c1 = 0
    cols = ["gene", "# sites", "# sites polymorphes", "# Syn", "# Non-Syn"]
    cols.extend(["# " + "-".join(cmp) for cmp in cmps])
    for col in cols:
        ws1.write(l1, c1, col, fmt_title); c1 += 1
    l1 += 1; c1 = 0

    # write titles to sheet2
    l2 = 0; c2 = 0
    cols = ["gene", "pos", "NT", "AA", "type"]
    cols.extend(group_list)
    cols.append("# selections")
    cols.extend(["-".join(cmp) for cmp in cmps])
    for col in cols:
        ws2.write(l2, c2, col, fmt_title); c2 += 1
    l2 += 1; c2 = 0

    # initialise summary counters
    summary = dict()
    for i in ['genes', 'sites', 'polymorphes_sites', 'S', 'NS']:
        summary[i] = 0
    summary['s_sites'] = dict()
    summary['s_genes'] = dict()
    for cmp1 in cmps:
        summary['s_sites'][cmp1] = 0
        summary['s_genes'][cmp1] = 0

    # for each fasta file (ie gene)
    for gene_file in args.fasta:
        gene_name = gene_file.split('/')[-1].replace('_aln2_NT.fasta', '').replace('_other.fasta', '')

        # detect variants
        sheet1_cols, sheet2_lines = detect_variants(gene_name, gene_file, groups, group_list, cmps, summary)

        # write count summary to sheet1
        for col in sheet1_cols:
            ws1.write(l1, c1, col); c1 += 1
        l1 += 1; c1 = 0

        # write sites polymorphic sites details to sheet2
        for line in sheet2_lines:
            for col in line:
                ws2.write(l2, c2, col); c2 += 1
            l2 += 1; c2 = 0

    # write per site summary to sheet1
    ws1.write(l1, c1, "sum sites", fmt_title); c1 += 1
    for i in ['sites', 'polymorphes_sites', 'S', 'NS']:
        ws1.write(l1, c1, summary[i], fmt_title); c1 += 1
    for cmp1 in cmps:
        ws1.write(l1, c1, summary['s_sites'][cmp1], fmt_title); c1 += 1
    l1 += 1; c1 = 0

    # write per gene summary to sheet1
    ws1.write(l1, c1, "sum {} genes".format(summary['genes']), fmt_title); c1 += 1
    for i in ['sites', 'polymorphes_sites', 'S', 'NS']:
        c1 += 1
    for cmp1 in cmps:
        ws1.write(l1, c1, summary['s_genes'][cmp1], fmt_title); c1 += 1

    # write excel file
    wb.save(args.xls)

def detect_variants(gene_name, gene_file, groups, group_list, cmps, summary):

    # read alignment file
    size, m_seq = read_aln(gene_name, gene_file)
    
    # select polymotphic sites and return matrix for codons and AA
    m_NT = select_polymorphes(m_seq, size)

    # detect non synonymous
    annot = annot_sites(gene_name, m_NT, groups, gene_file)

    # selected sites 
    # select[pos][cmp] <= 0/1
    # select[pos]['sum'] <= number of cmp selected for this pos
    # select['tot_selected']['cmp'] <= total number of variants selected from a comparison
    select = dict()
    select['tot_selected'] = dict()
    for cmp1 in cmps:
        select['tot_selected'][cmp1] = 0

    for pos in annot:
        select[pos] = dict()
        select[pos]['sum'] = 0
        for cmp1 in cmps:
            select[pos][cmp1] = "0"
    for cmp1 in cmps:
        for pos in select_variants(cmp1, groups, annot, m_NT):
            select[pos][cmp1] = "1"
            select[pos]['sum'] += 1
            select['tot_selected'][cmp1] += 1

    # summary of genes, sites, polymorphes sites, synonymous and non synonymous sites
    summary['genes'] += 1
    summary['sites'] += size
    summary['polymorphes_sites'] += len(m_NT)
    summary['S'] += annot['S']
    summary['NS'] += annot['NS']
    # per comparison summary of selected variants and genes
    for cmp1 in cmps:
        if select['tot_selected'][cmp1] > 0:
            summary['s_sites'][cmp1] += select['tot_selected'][cmp1] 
            summary['s_genes'][cmp1] += 1

    # summary count for sheet1 
    sheet1_cols = [gene_name, size, len(m_NT), annot['S'], annot['NS']]
    for cmp1 in cmps:
        sheet1_cols.append(select['tot_selected'][cmp1])
        
    # sites details for sheet2 
    sheet2_lines = list()
    for pos in sorted(annot):
        if not isinstance(pos, int):
            continue
        cols = [gene_name, pos, ",".join(annot[pos]['NT_alleles']), ",".join(annot[pos]['AA_alleles']), annot[pos]['type']]
        # pop details
        for pop in group_list:
            cols.append(",".join(annot[pos]['pops'][pop]))
        # cmp details
        cols.append(select[pos]['sum'])
        if select[pos]['sum'] > 0:
            for cmp1 in cmps:
                cols.append(select[pos][cmp1])
        sheet2_lines.append(cols)

    return sheet1_cols, sheet2_lines
        
def select_variants(cmp1, groups, annot, m_NT):
    """
    Select variants fixed in pop1 and pop2 and diff betwenn pop1 and pop2
    Returns a list of selected sites pos
    """

    pop1, pop2 = cmp1
    select = list()
    for pos in m_NT:
        if len(annot[pos]['pops'][pop1]) != 1:
            # site not fixed for pop1
            continue
        pop1_allele = annot[pos]['pops'][pop1][0]

        if len(annot[pos]['pops'][pop2]) != 1:
            # site not fixed for pop2
            continue
        pop2_allele = annot[pos]['pops'][pop2][0]

        if pop1_allele == pop2_allele:
            # site with same allele value in both pops
            continue

        # this site is selected : fixed in both pop and different betwenn pop
        # print "pos {} : diff between pops {} [{}] and {} [{}]".format(pos, pop1, pop1_allele, pop2, pop2_allele)
        select.append(pos)
    return select

def annot_sites(gene_name, m_NT, groups, gene_file):
    """
    annotations for each site (ie by position)
    1) per pop : alleles
       m_annot[pos]['pops'][pop] <= alleles list
    2) global : NT_alleles, AA_alleles, S/NS
       m_annot[pos]['NT_alleles'] <= list
       m_annot[pos]['AA_alleles'] <= list
       m_annot[pos]['type'] <= "S" or "NS"
       m_annot['S'] <= number of synonymous sites
       m_annot['NS'] <= number of non-synonymous sites

    !! for tRNA (trn*) and rRNA (rrn*) genes, these annotations are not set
       m_annot[pos]['AA_alleles'] <= list
       m_annot[pos]['type'] <= "S" or "NS"
       m_annot['S'] <= number of synonymous sites
       m_annot['NS'] <= number of non-synonymous sites
    """

    m_annot = dict()
    m_annot['S'] = 0
    m_annot['NS'] = 0

    for pos in m_NT:
        m_annot[pos] = dict()

        # 1) per pop : alleles
        m_annot[pos]['pops'] = dict()
        for pop in groups:
            alleles = set()
            for sample in groups[pop]:
                alleles.add(m_NT[pos][sample])
            m_annot[pos]['pops'][pop] = sorted(list(alleles))

        # 2) global : NT_alleles, AA_alleles, S/NS
        m_annot[pos]['NT_alleles'] = m_annot[pos]['pops']['all']

        # AA annot not needed for tRNA and rRNA genes
        if gene_name.startswith('trn') or gene_name.startswith('rrn'):
            m_annot[pos]['AA_alleles'] = list()
            m_annot[pos]['type'] = ""
        else:
            AA_list = list()
            for NT in m_annot[pos]['NT_alleles']:
                if not '-' in NT and not 'N' in NT:
                    AA_list.append("{}".format(Seq(NT).translate(table=11)))
            m_annot[pos]['AA_alleles'] = AA_list
            if len(set(AA_list)) == 1:
                m_annot[pos]['type'] = "S"
                m_annot['S'] += 1
            else:
                m_annot[pos]['type'] = "NS"
                m_annot['NS'] += 1

    if args.verbose:
        print "{} S {} NS {}".format(m_annot['S'], m_annot['NS'], gene_file)
    return m_annot

def select_polymorphes(m_seq, size):
    """
    Select polymorphes sites from all samples
    Returns a matrix (pos x sample) for polymorphes sites
    m_NT[pos][sample] <= 3 nuclueotids codon allele at this position for this sample
    """

    # list of codons positions (start=0, end=size, step=3) => 3, 6, 9 ...
    positions = range(0, size, 3)
    samples = sorted(m_seq)

    m_NT = dict()
    for pos in positions:
        check = set()
        NT = dict()
	# get codons for all individuals
	for sample in samples:
            codon = m_seq[sample][pos:(pos+3)]
            codon_str = "{}".format(codon)
            NT[sample] = codon_str
            if '-' not in codon_str and 'N' not in codon_str:
                check.add(codon_str)
	# test for polymorphism
	if len(check) > 1:
            # add to codon matrix
            m_NT[pos] = NT

    if args.verbose:
        print "{} / {} plymorphes sites".format(len(m_NT), len(positions))
    return m_NT

def read_aln(gene_name, gene_file):
    """
    Read multi fasta file
    Check that first seq's size is a multiple of 3 (execpt for tRNA (trn*) and rRNA (rrn*) genes)
    Check all seq have the same size that the first one
    !!!! WARNING for test : add - to smaller seq
    Skip sequences S* (outgroups) so that we keep only Silene nutans samples)
    return a table of seq, + seq size
      aln[sample] = sequence
    """

    aln = dict()
    size = 0

    # record the alignement (fasta)
    infile = parse(gene_file, 'fasta')
    for record in infile:
        if size == 0:
            size = len(record)
            # skip tRNA and rRNA
            if gene_name.startswith('trn'):
                continue
            if gene_name.startswith('rrn'):
                continue
            if size % 3 != 0:
                print "ERROR {} seq {} size {} not a multiple of 3".format(gene_file, record.id, len(record))
                exit(1)
        else:
            if len(record) != size:
                if len(record) != size:
                    record.seq += "-" * (size - len(record))
                else:
                    print "ERROR {} seq {} size {} != aln size {}".format(gene_file, record.id, len(record), size)
                    exit(1)

        # skip outgroups (name = S*)
	if record.id[0] == 'S':
            continue
            
        # clean seq name SAMPLE_GENE => SAMPLE
        sample = record.id.split('_')[0]

        aln[sample] = record.seq
    infile.close()

    if args.verbose:
        print "{} sites, {} samples from {}".format(size, len(aln), gene_file)
    return size, aln

def read_samples(samples_file):
    """ 
    read group sample from sample_file
    returns samples, groups
    samples[sample] = group
    groups[group] = sample_list
    groups = E, W, W1, W2, W3, all
    """
    
    samples = dict()
    tmp_groups = collections.defaultdict(list)
    groups = dict()
    with open(samples_file) as f:
        for l in f:
            if l.startswith('Ind'):
                continue
            if l.startswith('#'):
                continue
            ind, kapa, P5, P7, group = l.strip('\n').split(',')
            sample = kapa
            samples[sample] = group
            tmp_groups[group].append(sample)
            if group.startswith('W'):
                tmp_groups['W'].append(sample)
            tmp_groups['all'].append(sample)
            
    for group in tmp_groups:
        groups[group] = sorted(tmp_groups[group])
    groups_details = ", ".join(["{}={}".format(group, len(groups[group])) for group in groups])
    if args.verbose:
        print "{} sample, {}".format(len(samples), groups_details)
    return samples, groups

if __name__ == '__main__':
    description = """ 
    Detect variants from multi-fasta alignemnt 
    Determine synonymous or non synonymous variants type (using NCBI translation table 11 for plastid)
    Select variants seggregation for each comparison between 2 pops
    Create xls file with 2 sheets for summary count and polymorphic sites details
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--verbose", action="store_const", const="-v", default="")
    parser.add_argument("--fasta", nargs="+", help="input fasta files")
    parser.add_argument("--samples", help="input population file", default="samples.csv")
    parser.add_argument("--xls", help="output xls file")

    args = parser.parse_args()
    main() 
