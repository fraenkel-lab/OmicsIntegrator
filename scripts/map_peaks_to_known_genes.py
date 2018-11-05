#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import sys, os
import new_avl as avl
import heapq
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
parser.add_option('--upstream-window',dest='upst_win',type='int',default= 2000,help='window width in base pairs to consider promoter region [default: %default]')
parser.add_option('--downstream-window',dest='dnst_win',type='int',default= 2000,help='window width in base pairs to consider downstream region [default: %default]')
parser.add_option('--tss',dest='tss',action='store_true',default=True, help='calculate downstream window from transcription start site instead of transcription end site')
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
    gene_ref = dd(list)             #all of the genes in a chromosome in dictionary with keys = chromID and value = [list of genes]
    gene_info = {}                  #all information about a each gene, keys = geneID and value = dictunary of all info about gene
    chrom_info = {}                 #all of the chromosomes with info about gene stored in AVL tree, i.e keys - chromID and value = AVL(all genes arranged by startSyt)
    for ref_dict in reader :
        if ext.lower()=='.bed':
            ref_dict['txStart']=ref_dict['chromStart']
            ref_dict['txEnd']=ref_dict['chromEnd']

# determine intervals for promoter, gene, and downstream
        if  ref_dict['strand'] == '+' : #if gene in 5' to 3' orientation
            promoter_coords = max(ref_dict['txStart']-1-opts.upst_win,0), ref_dict['txStart']-1 #find the start and end of the promoter
            gene_coords = ref_dict['txStart'], ref_dict['txEnd'] #find the start and end of gene
            #use these coordinates if we're trying to window around TSS
            window_coords = ref_dict['txStart']+1,ref_dict['txStart']+opts.dnst_win
            downstream_coords = ref_dict['txEnd']+1, ref_dict['txEnd']+1+opts.dnst_win
            ref_dict['promoter_coords'] = promoter_coords
            ref_dict['gene_coords'] = gene_coords
            ref_dict['window_coords'] = window_coords
            ref_dict['downstream_coords'] = downstream_coords
        else :
            promoter_coords = ref_dict['txEnd']+1, ref_dict['txEnd']+1+opts.upst_win # +1 because we're using 1 based indexing
            gene_coords =ref_dict['txStart'], ref_dict['txEnd']
            window_coords = ref_dict['txEnd']-opts.dnst_win,ref_dict['txEnd']
            downstream_coords = ref_dict['txStart']-1-opts.dnst_win, ref_dict['txStart']-1 # -1 because we're using 1 based indexing
            ref_dict['promoter_coords'] = promoter_coords
            ref_dict['gene_coords'] = gene_coords
            ref_dict['window_coords'] = window_coords
            ref_dict['downstream_coords'] = downstream_coords

        gene_ref[ref_dict['chrom']].append(ref_dict)
        gene_info[ref_dict['name']] = ref_dict
        #putting relevant information about the gene into our AVL tree based on the chromosome
        if ref_dict['chrom'] not in chrom_info.keys():
            chrom_info[ref_dict['chrom']] = avl.AVLTree() #making a new instance of an AVL tree
            if ref_dict['strand'] == '+':
                chrom_info[ref_dict['chrom']].insert((ref_dict['promoter_coords'][0],ref_dict['downstream_coords'][1],ref_dict['name']))
            else :
                chrom_info[ref_dict['chrom']].insert((ref_dict['downstream_coords'][0],ref_dict['promoter_coords'][1],ref_dict['name']))
        else :
            if ref_dict['strand'] == '+':
                chrom_info[ref_dict['chrom']].insert((ref_dict['promoter_coords'][0],ref_dict['downstream_coords'][1],ref_dict['name']))
            else :
                chrom_info[ref_dict['chrom']].insert((ref_dict['downstream_coords'][0],ref_dict['promoter_coords'][1],ref_dict['name']))

    return gene_ref, gene_info, chrom_info

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

    gene_ref, gene_info, chrom_info = parse_gene_ref(args[0])
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
    totalrows=len(list(peaks_reader))
    peaks_reader = peaks_reader_cls(peaks_fn)#reopen to iterate
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


    peak_info = {} #all of the information about a peak in a dictionary. keys = peak_name and value = dictionary of all valuable information
    chrom_peaks = {} #all of the peaks in a chromosome in dictionary with keys = chromID and value = [list of genes]
    peak_number = 0
    for peak in peaks_reader:
        peak_number += 1
        # if this is a comment or header line get skip it
        #removed 'startswith' call so that this can work with tuples
        if peak[fieldnames[0]][0]=='#' or \
           peak[fieldnames[0]] == fieldnames[0] or \
           peak[fieldnames[0]][0]=='track' : continue

        for k,v in peak.items() : peak[k] = parse_number(v)

        # MACS output gives us summit
        if opts.peaks_fmt == 'MACS' :
            peak_loc = int(peak[start_field])+int(peak['summit'])
            peak['peakLoc'] = peak_loc
        elif opts.peaks_fmt == 'GPS' :
            #get position and also add in a real window
            ch,mid,tot = peak['Position']##reader already parses this into tuple
            peak['Position']=tot
            peak_loc = int(mid)
            peak['peakLoc'] = peak_loc
            peak[chr_field] = ch
            peak[start_field] = peak_loc-125
            peak[end_field] = peak_loc+125

        else : # peak assumed to be in the middle of the reported peak range
            peak_loc = (int(peak[start_field])+int(peak[end_field]))/2
            peak['peakLoc'] = peak_loc

        if peak['name'] == '.' or type(peak['name']) is int:     #if the peak does not have a name then assign it one based on the number of peaks we have seen so far
            peak['name'] = 'Peak_number_' + str(peak_number)
            if peak_number % 10000 == 0:
                print peak['name']

        if peak['name'] in peak_info.keys():            #if two peaks have the same ID, change it so they are not identical
            peak['name'] = (peak['name']) + 'd'

        peak_info[peak['name']] = peak
        if peak['chrom'] not in chrom_peaks.keys():
            chrom_peaks[peak['chrom']] = [(peak[start_field], peak['name'])]
        else:
            chrom_peaks[peak['chrom']].append((peak[start_field], peak['name']))

    peaks_writer = DictWriter(peak_output,output_fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
    peaks_writer.writerow(dict([(k,k) for k in output_fields]))
    unique_genes = set()
    map_stats = dd(int)
    rowcount=0

    interval=1000
    if totalrows > 100000:
        interval = 10000

    print '\nParsing %d rows from peak file and will provide update every %d rows'%(totalrows,interval)

    peaks_without_genes = []
    genes_without_peaks = []
    #walk through the peaks in a chromosome
    for chrom in chrom_peaks:
        heapq.heapify(chrom_peaks[chrom]) #sort them based on order on the chromosome
        while chrom_peaks[chrom] != []:
            pk = heapq.heappop(chrom_peaks[chrom]) #get the earliest appearing peak, thusfar, and remove it
            peak = peak_info[pk[1]]                #this is now the peak we will be looking at

            rowcount+=1
            if rowcount % interval ==0:
                print 'Processing row %d out of %d...'%(rowcount,totalrows)

            if chrom not in chrom_info:
                sys.stderr.write('WARNING: peak chromosome %s not found in gene reference, skipping: %s\n'%(peak[chr_field],peak))
                continue

            mapped = False

            #find all of the genes within the peak region
            genes_in_window = []
            first_gene = chrom_info[chrom].get_smallest_at_least(peak[start_field]) #find the first gene that begins after the start site of the peak
            if first_gene != None:
                node = first_gene
            else:
                node = chrom_info[chrom].find_biggest(chrom_info[chrom].rootNode)

            count = 100 #arbitrary large number
            #check the previous genes(genes that start before the start site of the peak) and see if they end after the start site
            while node is not None:
                if type(node.value) is list:    #this means there are multiple that start on the same start site
                    for gen in node.value:
                        if gen[1] >= peak[start_field]:
                            genes_in_window.append(gen)
                    node = chrom_info[chrom].predecessor(node) #set node to be the previous node

                elif node.value[1] < peak[start_field]: #if the end is after the start site of the peak, decrease count and go on till count = 0, ie there are 100 genes that end before the start of the peak
                    count -= 1
                    node = chrom_info[chrom].predecessor(node)
                    if count == 0:
                        break
                else:
                    genes_in_window.append(node.value)
                    node = chrom_info[chrom].predecessor(node)

            genes_in_window.reverse() # so they are in increasing not decreasing order

            #we should have gotten every gene that starts before a peak start site but ends before the peak ends
            #now we are going to look at all of the genes that begin after the start site of the peak till the end of the peak
            for elem in chrom_info[chrom].inorder(first_gene, peak[end_field]):
                if elem in genes_in_window:
                    continue
                elif type(elem) is list:
                    for i in elem:
                        genes_in_window.append(i)
                else:
                    genes_in_window.append(elem)

            if genes_in_window == []:   #we have found no genes in the peak region
                #sys.stderr.write('No genes in peak region, skipping: %s\n' %(peak['name']))
                peaks_without_genes.append(peak['name'])
                continue

            #loop through all of the genes that we have found to be relevant and mark them
            for genes in genes_in_window:
                # reusable dictionary for output
                out_d = {}.fromkeys(output_fields,0)
                out_d.update (peak)
                out_d['map type'] = ''
                out_d['chromo'] = peak[chr_field]
                if peak['peakLoc'] is None:
                    continue
                out_d['peak loc'] = peak['peakLoc']
                peak_loc = peak['peakLoc']

                gene = gene_info[genes[2]]
                promoter_coords = gene['promoter_coords']
                gene_coords = gene['gene_coords']
                window_coords = gene['window_coords']
                downstream_coords = gene['downstream_coords']

                # check for promoter
                if peak_loc >= promoter_coords[0] and peak_loc <= promoter_coords[1] :                      #if the peak location is in the promoter region
                    out_d['map type'] = 'promoter'
                    out_d['dist from feature'] = peak_loc - promoter_coords[1] if gene['strand'] == '+' else promoter_coords[0] - peak_loc

                ##now check to see if we are ONLY looking for window around tss
                # if gene['name'] == 'uc021wvh.1':
                #     print 'wtf'

                elif opts.tss:                                                                              #check if peak location is in transcription start site
                    if peak_loc >= window_coords[0] and peak_loc <= window_coords[1]:
                        out_d['map type'] = 'near TSS'
                        out_d['dist from feature'] = peak_loc - window_coords[0] if gene['strand'] == '+' else window_coords[1]-peak_loc
                    else:
                        continue
                # now we check for upstream/downstream areas if we are not using the opts.tss
                elif peak_loc >= gene_coords[0] and peak_loc <= gene_coords[1] :                            #check if peak location is within gene coordinates
                    # check for intron/exon
                    exon_coords = zip(gene['exonStarts'],gene['exonEnds'])                                  #make master list of tuples (exonStart, exonend)
                    in_exon = False
                    for st,en in exon_coords :
                        if peak_loc >= st and peak_loc <= en :                                              # check if the peak location is in an exon region
                            in_exon = True
                            break                                                                           #if peak location in exon, then exit out of loop
                    out_d['map type'] = 'gene'                                                              #map to gene anyways
                    out_d['map subtype'] = 'exon' if in_exon else 'intron'                                  #mark as exon if in_exon == true

                    # score = (peak-TSS)/(TSE-TSS) - peak distance from TSS as fraction of length of gene
                    gene_len = float(gene_coords[1]-gene_coords[0])                                         #find total length of gene
                    out_d['gene pos'] = (peak_loc-gene_coords[0])/gene_len if gene['strand'] == '+' else (gene_coords[1]-peak_loc)/gene_len #calculate peak distance from TSS as fraction of length of gene

                    # distance calculated from start of gene                                                #find the peak location from the end of the promoter end
                    out_d['dist from feature'] = peak_loc - promoter_coords[1] if gene['strand'] == '+' else promoter_coords[0] - peak_loc

                    map_stats[out_d['map subtype']] += 1

                # check for downstream if we're not doing a window
                elif peak_loc >= downstream_coords[0] and peak_loc <= downstream_coords[1] :                #if peak location is within downstream coordinates
                    out_d['map type'] = 'after'
                    out_d['dist from feature'] = peak_loc - downstream_coords[0] if gene['strand'] == '+' else downstream_coords[1] - peak_loc      #caluclate how far peak location is from downstream_coords

                # does not map to this gene
                else :                                                                                      #if you cannot find it within this gene
                    genes_without_peaks.append((gene['name'], peak['name']))
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

    if peak_output != sys.stdout:
        peak_output.close()


    #if opts.stats_output != sys.stderr :
    #    opts.stats_output = open(opts.stats_output,'w')

    #for k,v in map_stats.items() :
    #    opts.stats_output.write('%s: %s\n'%(k,v))

    #if opts.stats_output != sys.stderr :
    #    opts.stats_output.close()
