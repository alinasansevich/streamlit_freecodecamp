#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:16:09 2021

@author: alina

GC Content Calculator Web App
"""
import numpy as np
import pandas as pd
import streamlit as st
import altair as alt
from PIL import Image

# Page Title CONTINUE LATER XXXXXXXXXXXXXXXXXXXX
image = Image.open('GC-content_calculator.jpeg')

st.image(image, use_column_width=True)

st.write("""
# GC Content Calculator Web App

This app calculates the GC content of query DNA! XXXXXXX

***
""")

# Input Text Box

#st.sidebar.header('Enter DNA sequence')
st.header('Enter DNA sequence')
st.subheader('Paste sequence in FASTA format.')
st.write('Base Type: Only accept four letters ATGC (case-insensitive)')
st.write('Window Size: 30 nucleotides')

sequence_input = """>NM_119948.4 Arabidopsis thaliana phosphoenolpyruvate carboxykinase 1 (PCK1), mRNA
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
TACATGATAGAACCAAGTACGAATAATATAAATTCTGGTCCAATCTGATGATGATTTTCAAAA
"""

# sequence = st.sidebar.text_area("Sequence input", sequence_input, height=250)
sequence = st.text_area("Sequence input", sequence_input, height=250)
sequence = sequence.splitlines()
sequence = sequence[1:] # Skips the sequence name (first line)
sequence = ''.join(sequence) # Concatenates list to string

st.write("""
***
""")

#####################################
# What is GC Content?

# GC content is usually calculated as a percentage value and sometimes called G+C ratio or GC-ratio. GC-content percentage is calculated as Count(G + C)/Count(A + T + G + C) * 100%.

# Why to care about GC Content?

# The GC pair is bound by three hydrogen bonds, while AT pairs are bound by two hydrogen bonds. And so,
# The GC content affects the stability of DNA.
# The GC content affects the secondary structure of mRNA.
# The GC content affects the annealing temperature for template DNA in PCR experiments.
# Where to apply the GC Content?

# The GC Content can be used in
# 1. primer design for PCR experiments
# 2. gene design from protein expression
# 3. mRNA hairpin prediction
# and etc.

#####################################
# Count(G + C)/Count(A + T + G + C) * 100%
def gc_content(window):
    """
    (str) -> int
    Returns the GC content of window (a DNA or RNA sequence).
    """
    gc = window.count('G') + window.count('C')
    atgc = window.count('A') + window.count('T') + window.count('G') + window.count('C')
    
    return (gc/atgc) * 100

# "Walk" over the sequence and calculate the GC content for each 30 nucleotide window
all_gc = []
all_windows = []
for i in range(len(sequence) - 30):
    seq = sequence[i:i+30]
    gc_cont = gc_content(seq)
    all_gc.append(gc_cont)
    all_windows.append(seq)

# create GC content distribution graph
window_num = np.arange(1, len(all_gc)+1)
df_gc = pd.DataFrame({
    'window_num': window_num,
    'gc': all_gc
    })
gc_plot = alt.Chart(df_gc).mark_line(point=True).encode(x='window_num', y='gc')   
st.write(gc_plot)



# Display GC Content vs Window number
st.subheader('Window number vs GC content')
# df = pd.DataFrame.from_dict(X, orient='index')
# df = df.rename({0: 'count'}, axis='columns')
# df.reset_index(inplace=True)
# df = df.rename(columns = {'index':'nucleotide'})

################## round(x, 2)

# h_letters = [ letter for letter in 'human' ]

# all_gc_round = [round(x, 2) for x in all_gc] 

all_gc_round = []
for num in all_gc:
    all_gc_round.append(round(num, 2)

df_gc_windows = pd.DataFrame({
    'window_num': window_num,
    'gc': all_gc_round,
    'window': all_windows
    })

st.write(df_gc_windows)

### 1. Print dictionary
st.subheader('1. Print dictionary')
def DNA_nucleotide_count(seq):
  d = dict([
            ('A',seq.count('A')),
            ('T',seq.count('T')),
            ('G',seq.count('G')),
            ('C',seq.count('C'))
            ])
  return d

X = DNA_nucleotide_count(sequence)
### 3. Display DataFrame
st.subheader('3. Display DataFrame')
df = pd.DataFrame.from_dict(X, orient='index')
df = df.rename({0: 'count'}, axis='columns')
df.reset_index(inplace=True)
df = df.rename(columns = {'index':'nucleotide'})
st.write(df)
### 4. Display Bar Chart using Altair
st.subheader('4. Display Bar chart')
p = alt.Chart(df).mark_bar().encode(
    x='nucleotide',
    y='count'
)
p = p.properties(
    width=alt.Step(80)  # controls width of bar.
)
st.write(p)



# add something like this:
# Summary: Full Length(2583bp) | A(28% 724) | T(29% 739) | G(23% 598) | C(20% 522)
st.write("""
Summary: Full Length(2583bp) | A(28% 724) | T(29% 739) | G(23% 598) | C(20% 522)
***
""")


# gc content distribution graph
