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

pdb_paths=['/home/dell/Desktop/CPDAS/venv/pdb-files/P0A6P5_4_119.pdb', '/home/dell/Desktop/CPDAS/venv/pdb-files/P0A6P5_204_322.pdb', '/home/dell/Desktop/CPDAS/venv/pdb-files/P50743_5_120.pdb', '/home/dell/Desktop/CPDAS/venv/pdb-files/P50743_177_295.pdb']
alignment = AlignIO.read("domain_seq_aligned.aln", "clustal")
ids = []
for a in alignment.__dict__['_records']:
    ids.append(a.id)
id_path_dict = dict(zip(ids, pdb_paths))
PATH = "/home/dell/WebDrivers/chromedriver"
chrome_options = webdriver.ChromeOptions()
preferences = {
    "profile.default_content_settings.popups": 0,
    "download.default_directory": os.getcwd() + os.path.sep,
    "directory_upgrade": True
}
s = Service(PATH)
driver = webdriver.Chrome(service=s, options=chrome_options)
msa_path = str(os.getcwd()+r'/domain_seq_aligned.aln')
for key, value in id_path_dict.items():
        seq_id = key
        path = value
        driver.implicitly_wait(75)
        driver.get("https://consurf.tau.ac.il/")
        aa_button = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[2]/label/input").click()
        adv_options = driver.find_element(by=By.XPATH, value='/html/body/div[3]/div/form/table/tbody/tr/td[2]/a').click()
        yes_structure = driver.find_element(by=By.XPATH, value='/html/body/div[3]/div/form/div[1]/div[1]/label/input').click()
        upload_structure = driver.find_element(by=By.NAME, value='pdb_FILE')
        upload_structure.send_keys(path)
        chain_identifier = driver.find_element(by=By.XPATH, value='/html/body/div[3]/div/form/div[1]/div[2]/input[4]').send_keys('A')
        MSA_yes = driver.find_element(by=By.XPATH, value='/html/body/div[3]/div/form/div[3]/div[1]/div[1]/label/input').click()
        upload_msa = driver.find_element(by=By.XPATH, value='/html/body/div[3]/div/form/div[3]/div[2]/input[1]').send_keys(msa_path)
        try:
            ddmenu = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[3]/div[2]/div[2]/select")
            ddmenu = Select(ddmenu)
            ddmenu.select_by_visible_text(str(seq_id))
        except NoSuchElementException:
            seq_text = driver.find_element(by=By.NAME, value="msa_SEQNAME")
            seq_text.send_keys(str(seq_id))
        try:
            update_sel = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[3]/div[2]/div[2]/a[2]").click()
        except NoSuchElementException:
            link_text = driver.find_element(by=By.LINK_TEXT, value="update selection").click()
        try:
            no_email = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[13]/label/input").click()
        except NoSuchElementException:
            no_mail_check = driver.find_element(by=By.NAME, value="send_user_mail").click()
        try:
            submit_button = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/input[3]").click()
        except NoSuchElementException:
            driver.find_element(by=By.NAME, value="submitForm").click()
        try:
            alert = driver.switch_to.alert
            alert_msg = driver.switch_to.alert.text
            print(alert_msg)
            alert.accept()
        except NoAlertPresentException:
            print("exception handled")
        WebDriverWait(driver, 175).until(EC.presence_of_element_located((By.CLASS_NAME, "output_title")))
        WebDriverWait(driver, 20).until(EC.element_to_be_clickable((By.LINK_TEXT, "click!"))).click()
