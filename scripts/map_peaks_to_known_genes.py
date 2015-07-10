#!/usr/local/bin/python

import sys, os
from optparse import OptionParser
from collections import defaultdict as dd
from csv import DictReader, DictWriter



usage = '%prog [options] <knownGene file> <peaks file>'
description = """
Map the peaks in <peaks file> to genes in <knownGene file>.  <knownGene file> is\
format is as specified in http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.sql, though BED format is also accepted.\
<peaks file> format is as produced by GPS, MACS or BED.  If *auto* is chosen (default) file extension \
is examined for *.xls* for default MACS format, *.txt* for GPS, or *.bed* for BED format.  
"""
epilog = ''
parser = OptionParser(usage=usage,description=description,epilog=epilog)#,formatter=MultiLineHelpFormatter())
parser.add_option('--upstream-window',dest='upst_win',type='int',default=100000,help='window width in base pairs to consider promoter region [default: %default]')
parser.add_option('--downstream-window',dest='dnst_win',type='int',default=0,help='window width in base pairs to consider downstream region [default: %default]')
parser.add_option('--tss',dest='tss',action='store_true',default=False, help='calculate downstream window from transcription start site instead of transcription end site')
parser.add_option('--map-output',dest='peak_output',default=None,help='filename to output mapped peaks to [default: stdout]')
parser.add_option('--stats-output',dest='stats_output',default=sys.stderr,help='filename to output summary stats in conversion [default: stderr]')
parser.add_option('--peaks-format',dest='peaks_fmt',default='auto',type='choice',choices=['auto','MACS','BED'],help='format of peaks input file [default: %default]')
##TODO: add bigBed and GPS output

parser.add_option('--intergenic',dest='intergenic',action='store_true',help='write intergenic peaks to the gene file as well with None as gene ID')
parser.add_option('--utilpath',default='../src/',dest='addpath',help='Destination of chipsequtil library')
parser.add_option('--symbol-xref',dest='symbol_xref',default=None,help='Provide kgXref table file supplied to find a gene symbol and add as second column of output')

# TODO - options
#parser.add_option('--use-cds',dest='use_cds',action='store_true',help='use cdsStart and cdsEnd fields instead of txStart and txEnd to do mapping')
#parser.add_option('--capture-intergenic'...)
#parser.add_option('--map-format',dest='peak_format',type='choice',choices=['default','BED'],help='format of peak output [default: %default]')
#parser.add_option('--stats-format',dest='stats_format',type='choice',choices=['human','python'],help='format of summary stats output [default: %default]')

def parse_gene_ref(ref_gene) :
    ##some gene ids are not in txt form, provide ability to parse BED
    path,ext=os.path.splitext(ref_gene)
    if ext.lower()=='.bed':
        reader=BEDFile(ref_gene)
    else:
        reader = KnownGeneFile(ref_gene)
    gene_ref = dd(list)
    for ref_dict in reader :
        if ext.lower()=='.bed':
            ref_dict['txStart']=ref_dict['chromStart']
            ref_dict['txEnd']=ref_dict['chromEnd']
            
        gene_ref[ref_dict['chrom']].append(ref_dict)

    return gene_ref

def parse_gene_ref_line(l) :
    l = map(parse_number, l) # coerce to numbers where possible
    l[9] = map(parse_number, l[9].split(',')) # turn 'x,x,x,...' into list
    l[10] = map(parse_number, l[10].split(','))
    return l

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    sys.path.insert(0,opts.addpath)
    sys.path.insert(0,opts.addpath+'chipsequtil')
    
    from chipsequtil import MACSFile, BEDFile, KnownGeneFile, parse_number, GPSFile
 #  from chipsequtil.util import MultiLineHelpFormatter


    if len(args) < 2 :
        parser.error('Must provide two filename arguments')

    gene_ref = parse_gene_ref(args[0])
#    xref_fn = args[1]
    peaks_fn = args[1]
    if opts.peaks_fmt == 'auto' :
        path,ext = os.path.splitext(peaks_fn)
        if ext.lower() == '.xls' :
            opts.peaks_fmt = 'MACS'
        elif ext.lower() == '.bed' :
            opts.peaks_fmt = 'BED'
        elif ext.lower() == '.txt' :
            opts.peaks_fmt = 'GPS'
        else :
            parser.error('Could not guess peaks file format by extension (%s), aborting'%ext)

    if opts.peaks_fmt == 'MACS' :
        peaks_reader_cls = MACSFile
        chr_field, start_field, end_field = 'chr', 'start', 'end'
    elif opts.peaks_fmt == 'BED' :
        peaks_reader_cls = BEDFile
        chr_field, start_field, end_field = 'chrom', 'chromStart', 'chromEnd'
    elif opts.peaks_fmt == 'GPS' :
        peaks_reader_cls = GPSFile
        chr_field, start_field, end_field = 'chr', 'start', 'end'
    else :
        # should never happen
        fieldnames = []

    #peaks_reader = DictReader(open(args[1]),fieldnames=fieldnames,delimiter='\t')
    peaks_reader = peaks_reader_cls(peaks_fn)

    # default output format:
    if opts.peak_output :
        try:
            peak_output = open(opts.peak_output,'w')
        except:
            print "Error opening file:", sys.exc_info()[0]
            print "Check to make sure file exists at %s"%(opts.peak_output)
            raise
    else :
        peak_output = sys.stdout

    fieldnames = peaks_reader.FIELD_NAMES
    fieldnames += ["peak loc","dist from feature","gene pos","map type","map subtype"]
    output_fields = ['knownGeneID']+fieldnames

    # see if the user wants gene symbols too
    # TODO - actually make this an option, or make it required

    xref_fn=opts.symbol_xref
        
    if opts.symbol_xref :
        kgXref_fieldnames = ['kgID','mRNA','spID','spDisplayID','geneSymbol','refseq','protAcc','description']
        try:
            symbol_xref_reader = DictReader(open(opts.symbol_xref,'rU'),fieldnames=kgXref_fieldnames,delimiter='\t')
        except:
            print "Error opening file:", sys.exc_info()[0]
            print "Check to make sure file exists at %s"%(opts.symbol_xref)
            raise
        
        symbol_xref_map = {}
        for rec in symbol_xref_reader :
            symbol_xref_map[rec['kgID']] = rec
        output_fields = ['knownGeneID','geneSymbol']+fieldnames

    peaks_writer = DictWriter(peak_output,output_fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
    peaks_writer.writerow(dict([(k,k) for k in output_fields]))
    unique_genes = set()
    map_stats = dd(int)
    for peak in peaks_reader :

        # if this is a comment or header line get skip it
        #removed 'startswith' call so that this can work with tuples
        if peak[fieldnames[0]][0]=='#' or \
           peak[fieldnames[0]] == fieldnames[0] or \
           peak[fieldnames[0]][0]=='track' : continue

        # coerce values to numeric if possible
        for k,v in peak.items() : peak[k] = parse_number(v)

        # MACS output gives us summit
        if opts.peaks_fmt == 'MACS' :
            peak_loc = peak[start_field]+peak['summit']
        elif opts.peaks_fmt == 'GPS' :
            #get position and also add in a real window
            ch,mid,tot = peak['Position']##reader already parses this into tuple
            peak['Position']=tot
            peak_loc = int(mid)
            peak[chr_field] = ch
            peak[start_field] = peak_loc-125
            peak[end_field] = peak_loc+125
        else : # peak assumed to be in the middle of the reported peak range
            peak_loc = (peak[start_field]+peak[end_field])/2

        chrom_genes = gene_ref[peak[chr_field]]
        
        if len(chrom_genes) == 0 :
            sys.stderr.write('WARNING: peak chromosome %s not found in gene reference, skipping: %s\n'%(peak[chr_field],peak))
            continue

        mapped = False

        # walk through the genes for this chromosome
        for gene in chrom_genes :

            # reusable dictionary for output
            out_d = {}.fromkeys(output_fields,0)
            out_d.update(peak)
            out_d['map type'] = ''
            out_d['chromo'] = peak[chr_field]
            out_d['peak loc'] = peak_loc

            # determine intervals for promoter, gene, and downstream
            if gene['strand'] == '+' :
                promoter_coords = max(gene['txStart']-1-opts.upst_win,0), gene['txStart']-1
                gene_coords = gene['txStart'], gene['txEnd']
                #use these coordinates if we're trying to window around TSS
                window_coords = gene['txStart']+1,gene['txStart']+opts.dnst_win
               # if opts.tss :
                 #   gene_coords = gene['txStart'], min(gene['txEnd'],gene['txStart']+opts.dnst_win)
#                    downstream_coords = gene['txEnd']+1,gene['txStart']+opts.dnst_win
               #     downstream_coords = gene['txStart']+1,gene['txStart']+opts.dnst_win
               # else :
                #    gene_coords = gene['txStart'], gene['txEnd']
                downstream_coords = gene['txEnd']+1, gene['txEnd']+1+opts.dnst_win
            else :
                promoter_coords = gene['txEnd']+1, gene['txEnd']+1+opts.upst_win # +1 because we're using 1 based indexing
                gene_coords = gene['txStart'], gene['txEnd']
                window_coords = gene['txEnd']-opts.dnst_win,gene['txEnd'] 
                #if opts.tss :
                  #  gene_coords = max(gene['txStart'],gene['txEnd']-opts.upst_win), gene['txEnd']
#                    downstream_coords = gene['txEnd']-1-opts.dnst_win, gene['txStart']-1 # -1 because we're using 1 based indexing
                #    downstream_coords = gene['txEnd']-1-opts.dnst_win,gene['txEnd']-1
                #else :
                downstream_coords = gene['txStart']-1-opts.dnst_win, gene['txStart']-1 # -1 because we're using 1 based indexing

            # check for promoter
            if peak_loc >= promoter_coords[0] and peak_loc <= promoter_coords[1] :
                out_d['map type'] = 'promoter'
                out_d['dist from feature'] = peak_loc - promoter_coords[1] if gene['strand'] == '+' else promoter_coords[0] - peak_loc

            ##now check to see if we are ONLY looking for window around tss
            elif opts.tss:
                if peak_loc >= window_coords[0] and peak_loc <= window_coords[1]: 
                    out_d['map type'] = 'near TSS'
                    out_d['dist from feature'] = peak_loc - window_coords[0] if gene['strand'] == '+' else window_coords[1]-peak_loc
                else:
                    continue
            # now we can check for upstrema/downstream areas if we are not using the opts.tss
            elif peak_loc >= gene_coords[0] and peak_loc <= gene_coords[1] :
                # check for intron/exon
                exon_coords = zip(gene['exonStarts'],gene['exonEnds'])
                in_exon = False
                for st,en in exon_coords :
                    if peak_loc >= st and peak_loc <= en :
                        in_exon = True
                        break
                out_d['map type'] = 'gene'
                out_d['map subtype'] = 'exon' if in_exon else 'intron'

                # score = (peak-TSS)/(TSE-TSS) - peak distance from TSS as fraction of length of gene
                gene_len = float(gene_coords[1]-gene_coords[0])
                out_d['gene pos'] = (peak_loc-gene_coords[0])/gene_len if gene['strand'] == '+' else (gene_coords[1]-peak_loc)/gene_len

                # distance calculated from start of gene
                out_d['dist from feature'] = peak_loc - promoter_coords[1] if gene['strand'] == '+' else promoter_coords[0] - peak_loc

                map_stats[out_d['map subtype']] += 1

            # check for downstream if we're not doing a window
            elif peak_loc >= downstream_coords[0] and peak_loc <= downstream_coords[1] :
                out_d['map type'] = 'after'
                out_d['dist from feature'] = peak_loc - downstream_coords[0] if gene['strand'] == '+' else downstream_coords[1] - peak_loc

            # does not map to this gene
            else :
                pass

            # map type is not blank if we mapped to something
            if out_d['map type'] != '' :
                #out_d = {'knownGeneID':gene['name']}
                out_d['knownGeneID'] = gene['name']
                if opts.symbol_xref :
                    out_d['geneSymbol'] = symbol_xref_map[gene['name']]['geneSymbol']
                peaks_writer.writerow(out_d)

                mapped = True

                # reset map_type
                out_d['map type'] = ''

        if not mapped :
            if opts.intergenic :
                out_d['knownGeneID'] = 'None'
                out_d['geneSymbol'] = 'None'
                out_d['map type'] = 'intergenic'
                peaks_writer.writerow(out_d)
            map_stats['intergenic'] += 1

    if peak_output != sys.stdout :
        peak_output.close()

    #if opts.stats_output != sys.stderr :
    #    opts.stats_output = open(opts.stats_output,'w')

    #for k,v in map_stats.items() :
    #    opts.stats_output.write('%s: %s\n'%(k,v))

    #if opts.stats_output != sys.stderr :
    #    opts.stats_output.close()
