#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 15:10:59 2023
BOLD SCRAPER, but no SQL bullshit- I was making stuff way too complicated for the intended purposes...
Lets just write something that can easily be plugged in.
@author: christian
"""
import numpy as np
import subprocess as sbp
from Bio import AlignIO
from itertools import product
# base functions
# the actual meat
 #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# get dictionary of BOLD sequences 
def scrubBOLD(searchterm):
    import requests
    front = 'http://v3.boldsystems.org/index.php/API_Public/sequence?taxon='
    send = front + searchterm
    call = requests.get(send)
    call = call.content
    call = str(call)
    call = call.replace('\\r\\n', '\n')
    call = call.split('\n')
    call = dict(zip(call[::2], call[1::2]))
    return call

# concensus sequences
def conseq(aln):
    untransposed = list()
    concensus = list()
    for sequence in aln:
        untransposed.append(sequence.seq)
    transposed = np.array(untransposed).T.tolist()
    for i in transposed:
        max_nucleotide = max(set(i), key = i.count)
        if max_nucleotide != '-':
            concensus.append(max_nucleotide)
        elif max_nucleotide == '-':
            i_omitted_missing = [n for n in i if n != '-']
            if len(i_omitted_missing) == 0:
                concensus.append('-')
            elif len(i_omitted_missing) >= 1:
                max_nucleo_omit = max(set(i_omitted_missing), key = i_omitted_missing.count)
                concensus.append(max_nucleo_omit)
    con_seq = str()
    con_seq = con_seq.join(concensus)
    return con_seq
 #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
     
# =============================================================================
# =============================================================================
# some time of annotation
class bsequence(object):
    def __init__(self, key, sequence):
        self.key = key
        self.sequence = sequence
        self.boldID = self.reformatBOLDnames()[0]
        self.Genus_sp = self.reformatBOLDnames()[1]
        self.gene = self.reformatBOLDnames()[2]
        
    def reformatBOLDnames(self):
        boldname = self.key
        boldname = boldname.lstrip('>')
        boldname = boldname.split('|')
        return boldname
    
class csequence(object):
    def __init__(self, Genus_sp, gene, sequence):
        self.Genus_sp = Genus_sp
        self.gene = gene
        self.sequence = sequence
# =============================================================================
# =============================================================================
 # maybe some more annotation
class container(object):
    # __init__
    def __init__(self):
        self.sequences = list()
        self.current_subsection = self.sequences
        self.current_sp = str()
        self.current_gene = str()
        self.concensus_sequences = dict()
        self.genes = list()
        self.Genus_spp = list()
        self.temp_alignment_str = '_temp_align_.fa'
        self.current_fasta_filepath = str()
        self.temp_aligned_str = '_aligned.fa'
        self.current_aligned_path = str()
        self.current_alignment = 'no current alignments'
        self.__single_records__ = list()
    # -------------------------------------------------------------------------
    # add bsequences
    def add(self, bsequence):
        self.sequences.append(bsequence)
        
    def show(self):
        for bseq in self.current_subsection:
            print(f'{bseq.Genus_sp}, {bseq.gene}, {bseq.boldID}')

    # allgenes
    def allgenes(self):
        genes = list()
        for seq in self.sequences:
            if seq.gene in genes:
                pass
            else:
                genes.append(seq.gene)
        self.genes = genes
        return genes
    
    # allspecies
    def allspecies(self):
        species = list()
        for sp in self.sequences:
            if sp.Genus_sp in species:
                pass
            else:
                species.append(sp.Genus_sp)
        self.Genus_spp = species
        return species
    # -------------------------------------------------------------------------
    def subsection(self, species, gene):
        selection = list()
        for bseq in self.sequences:
            if bseq.Genus_sp == species and bseq.gene == gene:
                selection.append(bseq)
            else:
                pass
        self.current_subsection = selection
        self.current_sp = species
        self.current_gene = gene
    
    # writes current selection to fasta
    def fasta(self):
        filepath = f'{self.current_sp}_{self.current_gene}{self.temp_alignment_str}'.replace(' ',
                                                                                             '-')
        self.current_fasta_filepath = filepath
        for bseq in self.current_subsection:
            header = f'>{bseq.Genus_sp}_{bseq.gene}_{bseq.boldID}'
            body = bseq.sequence
            with open(filepath, 'a') as fasta:
                fasta.write(f'{header}\n{body}\n')
    # mafft -- # we should make the string passed to subprocess as an alternative
    # string with arguements which can be passed in order - easy to do, but we can do it later ...
    def mafft(self):
        aligned_path = f'{self.current_sp}_{self.current_gene}{self.temp_aligned_str}'.replace(' ',
                                                                                               '-')
        self.current_aligned_path = aligned_path
        mafft_command = '"/usr/bin/mafft"  --auto --clustalout --reorder {fa} > {out}'.format(fa = self.current_fasta_filepath,
                                                                                              out = self.current_aligned_path)
        sbp.call(mafft_command, shell = True)

    def openalignment(self, aln_format = 'clustal'):
        aln = AlignIO.read(self.current_aligned_path, aln_format)
        self.current_alignment = aln
        return aln
    
    # notice: this concensus is GREEDY; it prefers nucleotides to gaps.
    def concensus(self):
        untransposed = list()
        concensus = list()
        for sequence in self.current_alignment:
            untransposed.append(sequence.seq)
        transposed = np.array(untransposed).T.tolist()
        for i in transposed:
            max_nucleotide = max(set(i), key = i.count)
            if max_nucleotide != '-':
                concensus.append(max_nucleotide)
            elif max_nucleotide == '-':
                i_omitted_missing = [n for n in i if n != '-']
                if len(i_omitted_missing) == 0:
                    concensus.append('-')
                elif len(i_omitted_missing) >= 1:
                    max_nucleo_omit = max(set(i_omitted_missing), key = i_omitted_missing.count)
                    concensus.append(max_nucleo_omit)
        con_seq = str()
        con_seq = con_seq.join(concensus)
        self.concensus_sequences.update({f'{self.current_sp}_{self.current_gene}':con_seq})
        return csequence(self.current_sp, self.current_gene, con_seq)
    
    def concensus_combinations(self):
        species = self.allspecies()
        genes = self.allgenes()
        for sp_gene in product(species, genes):
            self.subsection(*sp_gene)
            self.fasta()
            self.mafft()
            try:
                self.openalignment()
                self.concensus()
            except ValueError:
                self.__single_records__.append(sp_gene)
                

# =============================================================================
# =============================================================================             

# DON'T RUN BELOW (use-case)
# %%
w = scrubBOLD('Chetogena scutellaris')
bigw = container()
for i,e in w.items():
    bigw.add(bsequence(i,e))
# %%
bigw.subsection('Chetogena parvipalpis', 'COI-5P')
bigw.fasta()
bigw.mafft()
bigw.openalignment()
bigw.concensus()
# %%
## OR
# %%
bigw.concensus_combinations()
# %%














