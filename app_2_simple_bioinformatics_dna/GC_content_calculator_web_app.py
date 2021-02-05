#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:16:09 2021

@author: alina

GC Content Calculator Web App

Idea from: https://www.youtube.com/watch?v=JwSS70SZdyM
"""
import numpy as np
import pandas as pd
import streamlit as st
import altair as alt
from PIL import Image


# Count(G + C)/Count(A + T + G + C) * 100%
def gc_content(window):
    """
    (str) -> int
    Returns the GC content of window (a DNA or RNA sequence).
    """
    gc = window.count('G') + window.count('C')
    atgc = window.count('A') + window.count('T') + window.count('G') + window.count('C')
    
    return (gc/atgc) * 100


def DNA_nucleotide_count(seq):
    """
    (str) -> dict
    Returns a dictionary with the nucleotide frequency for sequence seq.
    """
    d = dict([
        ('A',seq.count('A')),
        ('T',seq.count('T')),
        ('G',seq.count('G')),
        ('C',seq.count('C'))
        ])
    return d

def seq_summary(nucleotide_count):
    """
    (dict) -> str

    Returns a summary of theljkgjcgh
    """
    full_length = sum(nucleotide_count.values())
    a = round(nucleotide_count['A'] * 100 / full_length)
    t = round(nucleotide_count['T'] * 100 / full_length)
    g = round(nucleotide_count['G'] * 100 / full_length)
    c = round(nucleotide_count['C'] * 100 / full_length)
    
    summary = "Summary: Full Length({}bp) | A({}% {}) | T({}% {}) | G({}% {}) | C({}% {})".format(
        full_length, a, nucleotide_count['A'], t, nucleotide_count['T'], g, nucleotide_count['G'], c, nucleotide_count['C'],)
    
    return summary


##### Page Title + Text
image = Image.open('GC-content_calculator.jpeg')

st.image(image, use_column_width=True)

st.write("""
# GC Content Calculator Web App
***
What is GC Content?
 
GC content percentage is calculated as Count(G + C)/Count(A + T + G + C) * 100%.

The GC pair is bound by three hydrogen bonds, while AT pairs are bound by two hydrogen bonds.
This means that the GC content affects the stability of DNA molecules, the secondary structure of mRNAs
and the annealing temperature for primers and template DNA in PCR experiments.

Knowing the GC-content of a DNA region is useful when designing primers for PCR experiments,
since a higher GC-content level indicates a relatively higher melting temperature.
***
""")

##### Input Text Box
st.sidebar.header('GC Content Calculator Web App')    ### 
st.header('Enter DNA sequence')
st.subheader('Paste sequence in FASTA format, then press Ctrl+Enter to apply.')
st.write('Base Type: Only accepts four letters ATGC (case-insensitive)')
st.write('Window Size: 30 nucleotides')

sequence_input = """>Test sequence
AGTCATCTTTATAAACCACCGGTTATGTTAAGAGAGAAAATAAAAATAAAAAAGGGGCTCTTCCTAGGAA
GATAGATCTTAACCATGGTTAACACTCTCACGGTTCATTATTAACCATGGTTCTAAAAATCTAACCTTTA
AAAAACCACTTTCGCTTCTCTTCACATTCGCATCATTTTGTATCATCCCTTGAAAACGTTAAATGATCTT
CTCCTCCGATCATTAGTCTCTTTAATCTTTCTCAGCCTCTTCTTGTTCGTGATCTCTCTTCCTCCGGAAA
AAGATGTCGGCCGGTAACGGAAATGCTACTAACGGTGACGGAGGGTTTAGTTTCCCTAAAGGACCGGTGA
TGCCGAAGATAACGACCGGAGCAGCAAAGAGAGGTAGCGGAGTCTGCCACGACGATAGTGGTCCGACGGT
GAATGCCACAACCATCGATGAGCTTCATTCGTTACAGAAGAAACGTTCTGCTCCTACCACACCGATCAAC
CAAAACGCCGCCGCTGCTTTTGCCGCCGTCTCCGAGGAGGAGCGTCAGAAGATTCAGCTTCAATCTATCA
GTGCATCGTTAGCATCGTTAACGAGAGAGTCAGGACCAAAGGTGGTGAGAGGAGATCCGGCGGAGAAGAA
GACCGATGGTTCAACTACTCCGGCGTACGCTCACGGCCAACATCATTCTATCTTTTCTCCGGCTACTGGT
GCTGTCAGTGATAGCTCCTTGAAGTTTACTCACGTCCTCTACAATCTTTCGCCTGCAGAGCTTTATGAGC
AAGCTATTAAGTATGAGAAAGGTTCGTTTATCACTTCTAATGGAGCTTTGGCGACGCTTTCTGGTGCTAA
GACTGGTCGTGCTCCCAGAGATAAGCGTGTTGTTAGAGATGCTACTACTGAGGATGAGCTTTGGTGGGGA
AAGGGTTCGCCGAATATCGAAATGGATGAACATACTTTCATGGTGAACAGAGAAAGAGCTGTTGATTACT
TGAATTCCTTGGAAAAGGTCTTTGTCAATGACCAATACTTAAACTGGGATCCAGAGAACAGAATCAAAGT
CAGGATTGTCTCAGCTAGAGCTTACCATTCATTGTTTATGCACAACATGTGTATCCGACCAACTCAGGAG
GAGCTTGAGAGCTTTGGTACTCCGGATTTTACTATATACAATGCTGGGCAGTTTCCATGTAATCGTTACA
CTCATTACATGACTTCGTCCACTAGCGTAGACCTTAATCTGGCTAGGAGGGAAATGGTTATACTTGGTAC
TCAGTATGCTGGGGAAATGAAGAAGGGTCTTTTCAGTGTGATGCATTACCTTATGCCTAAGCGTCGTATT
CTCTCCCTTCATTCTGGATGCAATATGGGAAAAGATGGAGATGTTGCTCTCTTCTTTGGACTTTCAGGTA
CCGGGAAGACAACGCTGTCTACTGATCACAACAGGTATCTTATTGGAGATGATGAGCATTGTTGGACTGA
GACTGGTGTTTCGAACATTGAGGGTGGGTGCTATGCTAAGTGTGTTGATCTTTCGAGGGAGAAGGAGCCT
GATATCTGGAACGCTATCAAGTTTGGAACAGTTTTGGAAAATGTTGTGTTTGATGAGCACACCAGAGAAG
TGGATTACTCTGATAAATCTGTTACAGAGAACACACGTGCTGCCTACCCAATTGAGTTCATTCCAAATGC
GAAAATACCTTGTGTTGGTCCACACCCGACAAATGTGATACTTCTGGCTTGTGATGCCTTTGGTGTTCTC
CCACCTGTGAGCAAGCTGAATCTGGCACAAACCATGTACCACTTCATCAGTGGTTACACTGCTCTGGTTG
CTGGCACAGAGGATGGTATCAAGGAGCCAACAGCAACATTCTCAGCTTGCTTTGGTGCAGCTTTCATAAT
GTTGCATCCCACAAAGTATGCAGCTATGTTAGCTGAGAAGATGAAGTCACAAGGTGCTACTGGTTGGCTC
GTCAACACTGGTTGGTCTGGTGGCAGTTATGGTGTTGGAAACAGAATCAAGCTGGCATACACTAGAAAGA
TCATCGATGCAATCCATTCGGGCAGTCTCTTGAAGGCAAACTACAAGAAAACCGAAATCTTTGGATTTGA
AATCCCAACTGAGATCGAAGGGATACCTTCAGAGATCTTGGACCCCGTCAACTCCTGGTCTGATAAGAAG
GCACACAAAGATACTCTGGTGAAACTGGGAGGTCTGTTCAAGAAGAACTTCGAGGTTTTTGCTAACCATA
AGATTGGTGTGGATGGTAAGCTTACGGAGGAGATTCTCGCTGCTGGTCCTATCTTTTAGAAAAACCAAAC
TCTGAATGATGTTGTCGAAAAGAAAGAAAGATCTACTATTATTAAGAAGAATAAAATGAGACTTTGTGTT
TTTCTTTGCTGTGATAATCTCTCTGAATAATAGAGAGAATTTAATGTTCCAAATGTGGTGCGATATGGAT
AATGATGATGACTATATGTAATTTATTACTATCCATCGTCATATATTTTTGTTGGGCCTTTGCCAACATT
TACATGATAGAACCAAGTACGAATAATATAAATTCTGGTCCAATCTGATGATGATTTTCAAAA"""

# remove header
sequence = st.text_area("Sequence input", sequence_input, height=250)
sequence = sequence.splitlines()
sequence = sequence[1:] # Skips the sequence name (first line)
sequence = ''.join(sequence) # Concatenates list to string

st.write("""
***
""")

# "Walk" over the sequence and calculate the GC content for each 30 nucleotide window
all_gc = []
all_windows = []
for i in range(len(sequence) - 30):
    seq = sequence[i:i+30]
    gc_cont = gc_content(seq)
    all_gc.append(gc_cont)
    all_windows.append(seq)

##### create GC content distribution graph
window_num = np.arange(1, len(all_gc)+1)
df_gc = pd.DataFrame({
    'Window Number': window_num,
    'GC%': all_gc
    })

# plot
gc_plot = alt.Chart(df_gc).mark_line(point=True).encode(x='Window Number', y='GC%')   
st.write(gc_plot)

##### Display Summary:
nucleotide_count = DNA_nucleotide_count(sequence)
s = seq_summary(nucleotide_count)
st.write(s)
st.write("""***""")
# add something like this:
# Summary: Full Length(2583bp) | A(28% 724) | T(29% 739) | G(23% 598) | C(20% 522)


##### Display DataFrame: GC Content vs Window number
st.subheader('GC content for each window:')
# df = pd.DataFrame.from_dict(X, orient='index')
# df = df.rename({0: 'count'}, axis='columns')
# df.reset_index(inplace=True)
# df = df.rename(columns = {'index':'nucleotide'})

all_gc_round = [round(x) for x in all_gc] 

gc_windows = {'Window #': window_num,
              'GC content': all_gc_round,
              'Sequence': all_windows}
df_gc_windows = pd.DataFrame.from_dict(gc_windows)

st.write(df_gc_windows)
