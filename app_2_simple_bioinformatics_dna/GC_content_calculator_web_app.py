#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:16:09 2021

@author: alina

GC Content Calculator Web App
"""

import pandas as pd
import streamlit as st
import altair as alt
from PIL import Image

# Page Title
image = Image.open('GC-content_calculator.jpeg')

st.image(image, use_column_width=True)

st.write("""
# GC Content Calculator Web App

This app calculates the GC content of query DNA! XXXXXXX

***
""")

