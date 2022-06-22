import os.path
import subprocess
import subprocess
import pandas as pd
import time
from selenium import webdriver
from bs4 import BeautifulSoup
from Functions import *
from Bio import SeqIO
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
import zipfile
from six.moves import urllib
from pdb2sql import pdb2sql
import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib as plt
import seaborn as sns
import plotly.express as px
import ast

PDBNAME = input("Enter PDB name: ")
input_scores = input("Enter scores list: ")
range_scores = ast.literal_eval(input_scores)
#print(PDBNAME, range_scores, type(range_scores))

db1 = pdb2sql(PDBNAME)
#db1.print()
resseq = []
for x in db1.get('resSeq'):
    for i in x:
        resseq.append(i)
score_l = [range_scores[i] for i in resseq]
val = np.array(score_l, dtype=float)
db1.update_column('temp', values=val,)
#db1.print()
db1.exportpdb("{}_temp_modified.pdb".format(PDBNAME))
