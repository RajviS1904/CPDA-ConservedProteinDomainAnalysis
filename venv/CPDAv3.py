#Import Statements

from Functions import *
import os.path
import subprocess
import subprocess
import pandas as pd
import time
from bs4 import BeautifulSoup
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import zipfile
from six.moves import urllib
from pdb2sql import pdb2sql
import numpy as np
from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.support.select import Select
from selenium.webdriver.common.alert import Alert
from selenium.common.exceptions import NoAlertPresentException
from selenium.common.exceptions import NoSuchElementException


#Input the sequence FASTA file name. (Note that the sequence id should be in the following format: >sp|P0A6P5|DER_ECOLI)

name = str(input("Enter your file name: "))
fin = []
#Fetch domains information from PFAM:
Fetch_PFAM(fname=name)

#Slice sequences into respective protein domains
slicedseq_PFAMdata(fin)

#Enter familyID and sort out all the required protein domains from the input sequences
domain_dataframe("MMR_HSR1", df3)

#writeout the sliced sequences into a new file
sliceddomains_FASTA(df4)

#align domains
align_domains()

#format aligned files
format_alignout()

#calculate conservation scores
conservation_scores_calc()

#format conservation scores into dataframes
conservation_scores_df()

#format normalized scores into dataframes
normalized_scores_df()

#predict variable residues from the alignment
predicted_variable_res(df6)

#format variable residues into dataframes
