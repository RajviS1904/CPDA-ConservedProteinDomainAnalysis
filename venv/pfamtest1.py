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

fname = "testseq"
fin = []
for record in SeqIO.parse(fname+".fasta", "fasta"):
    rec_id = (">" + record.id)
    seq = record.seq
    PATH = "/home/dell/WebDrivers/chromedriver"
    chrome_options = webdriver.ChromeOptions()
    #chrome_options.add_argument('headless')
    #chrome_options.add_argument('window-size=1920x1080')
    #chrome_options.add_argument("disable-gpu")
    s = Service(PATH)
    driver = webdriver.Chrome(service=s, options=chrome_options)
    driver.implicitly_wait(75)
    driver.get("http://pfam.xfam.org/")
    rec_search = driver.find_element(by=By.XPATH, value='/html/body/div[5]/div[2]/ul/li[5]').click()
    recordid = record.id[3:9]
    searchbox = driver.find_element(by=By.XPATH, value='/html/body/div[5]/div[7]/form/div[1]/input[1]').send_keys(recordid)
    driver.find_element(by=By.XPATH, value='/html/body/div[5]/div[7]/form/div[2]/a[1]').click()
    timeout = 175
    try:
        element_present = EC.presence_of_element_located((By.ID, 'content'))
        WebDriverWait(driver, timeout).until(element_present)
    except TimeoutException:
        print("Timed out waiting for page to load")
        break
    html = driver.page_source
    soup = BeautifulSoup(html, "lxml")
    div = soup.find_all('table', {'id': 'imageKey'})
    table = pd.read_html(str(div))
    for t in table:
        t['id'] = rec_id
    fin.append(table)
    driver.close()
    print("Finished getting PFAM Results Table data")


def slicedseq_PFAMdata(fin):
    df3 = pd.DataFrame()
    for f in fin:
        df = f[0]
        df2 = df.T.reset_index(level=1).T
        df2.columns = [i + "-" + j if (i != j and j != "") else i for i, j in zip(df2.columns, df2.loc['level_1'])]
        df2 = df2.drop('level_1', axis=0)
        df3 = pd.concat([df3, df2])
    data = []
    for record in SeqIO.parse(fname + ".fasta", "fasta"):
        data.append([">" + record.id, "".join(record.seq)])
    df_seq = pd.DataFrame(data, columns=['id', 'seq'])
    df3 = df3.merge(df_seq, on='id')
    df3['sliced_sequence'] = df3.apply(lambda x: x['seq'][x['Envelope-Start'] - 1:x['Envelope-End']], axis=1)
    df3['new_id'] = df3.apply(lambda x: x['id'] + "|" + str(x['Envelope-Start']) + "-" + str(x['Envelope-End']), axis=1)
    print(df3)

slicedseq_PFAMdata(fin)
