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


#Input the sequence FASTA file name. (Note that the sequence id should be in the following format: >sp|P0A6P5|DER_ECOLI)
File_name = str(input("Enter your file name: "))
fin = []

for record in SeqIO.parse(File_name + ".fasta", "fasta"):
    rec_id = (">" + record.id)
    seq = record.seq
    PATH = "/home/dell/WebDrivers/chromedriver"
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('headless')
    chrome_options.add_argument('window-size=1920x1080')
    chrome_options.add_argument("disable-gpu")
    s = Service(PATH)
    driver = webdriver.Chrome(service=s, options=chrome_options)
    driver.implicitly_wait(75)
    driver.get("http://pfam.xfam.org/")
    seq_search = driver.find_element(by=By.XPATH, value='/html/body/div[5]/div[2]/ul/li[2]')
    seq_search.click()
    sequence = seq
    searchbox = driver.find_element(by=By.ID, value='seq').send_keys(sequence)
    driver.find_element(by=By.XPATH, value='/html/body/div[5]/div[4]/div/a[1]/span').click()
    timeout = 175
    try:
        element_present = EC.presence_of_element_located((By.ID, "pfamASummary"))
        WebDriverWait(driver, timeout).until(element_present)
    except TimeoutException:
        print("Timed out waiting for page to load")
        break
    html = driver.page_source
    soup = BeautifulSoup(html, "lxml")
    div = soup.find_all('table', {'id': 'pfamASummary'})
    table = pd.read_html(str(div))
    for t in table:
        t['id'] = rec_id
    fin.append(table)
    driver.close()
print("Finished getting PFAM Results Table data")

# Format PFAM Results using pandas Dataframe

df3 = pd.DataFrame()
for f in fin:
    df = f[0]
    df2 = df.T.reset_index(level=1).T
    df2.columns = [i + "-" + j if (i != j and j != "") else i for i, j in zip(df2.columns, df2.loc['level_1'])]
    df2 = df2.drop('level_1', axis=0)
    df3 = pd.concat([df3, df2])
data = []
for record in SeqIO.parse(File_name + ".fasta", "fasta"):
    data.append([">" + record.id, "".join(record.seq)])
df_seq = pd.DataFrame(data, columns=['id', 'seq'])
df3 = df3.merge(df_seq, on='id')
df3['sliced_sequence'] = df3.apply(lambda x: x['seq'][x['Envelope-Start'] - 1:x['Envelope-End']], axis=1)
df3['new_id'] = df3.apply(lambda x: x['id'] + "|" + str(x['Envelope-Start']) + "-" + str(x['Envelope-End']), axis=1)
ln = df3['Family'].count()
df4 = pd.DataFrame()
for x in range(ln):
    if df3['Family'].iloc[x].find("MMR_HSR1") == 0:
        df4 = pd.concat([df4, df3.iloc[x:x + 1, :]])
    else:
        continue
df4.reset_index(drop=True)

#Write sliced sequences to a new FASTA File
with open("sliced_sequences.fasta", "w") as file:
    for idx, row in df4.iterrows():
        file.write(row['new_id'] + "\n")
        file.write(row['sliced_sequence'] + "\n")
print("wrote sliced sequences in new file")

#Align sliced sequences using PROBCONS
align = "./probcons/probcons -clustalw sliced_sequences.fasta > sliced_seq_align.aln"
subprocess.call(align, shell=True)

#Update the first line of Clustal Output File:
with open('sliced_seq_align.aln', 'r') as file:
    data = file.readlines()
data[0] = "CLUSTAL W (1.82) multiple sequence alignment\n"
with open('sliced_seq_align.aln', 'w') as file:
    file.writelines(data)
#Read alignment file
alignment = AlignIO.read("sliced_seq_align.aln", "clustal")

l = []
for i in alignment.__dict__['_records']:
    l.append([i.id, "".join(i.seq)])
df5 = pd.DataFrame(l)
df5.columns = ['id', 'seq']
df6 = pd.DataFrame([list(word) for word in list(df5.seq)])
df7 = pd.DataFrame((df6[col].value_counts()) for col in range(len(alignment[0].seq)))
df8 = df7.agg(['idxmax', 'max'], axis=1).mask(lambda x: x['max'].eq(0))
df9 = df8.T
word = df9.values.tolist()
consensus = ''.join(word[0])
df10 = pd.DataFrame({'id': ['Consensus'], 'seq': [consensus]})
df11 = pd.concat([df5, df10], ignore_index=True, axis=0)
df11.set_index("id", inplace=True)
df12 = pd.DataFrame([list(word) for word in list(df11.seq)], index=df11.index)
list1 = [dict(zip([(df12.columns[x])], [(np.where(df12.iloc[x, :] == df12.loc["Consensus"], True, False)).tolist()]))
         for x in range(len(df11))]
d1 = {}
for i in list1:
    d1.update(i)
df13 = pd.DataFrame(d1)
df14 = df13.T
df14.index = df11.index
#Creating a true false matrix for sequence comparison with consensus
df14.replace({False: 0, True: 1}, inplace=True)

#Calculate conservation scores using AACons Java Executable
cons_score = "java -jar compbio-conservation-1.1.jar -i=sliced_seq_align.aln -o=cons_out.csv"
norm_score = "java -jar compbio-conservation-1.1.jar -i=sliced_seq_align.aln -n -o=norm_out.csv"
subprocess.call(norm_score, shell=True)
subprocess.call(cons_score, shell=True)

#Format conservation scores in dataframes
df15 = pd.read_csv('cons_out.csv', delimiter='\s+', header=None)
df15.set_index([0], inplace=True)
df15.reset_index
df15.index.rename('Method', inplace=True)
df15.columns = pd.RangeIndex(0, len(df15.columns))

#Format normalized conservation scores in dataframes
df16 = pd.read_csv('norm_out.csv', delimiter='\s+', header=None)
df16.set_index([0], inplace=True)
df16.reset_index
df16.index.rename('Method', inplace=True)
df16.columns = pd.RangeIndex(0, len(df16.columns))

NormScoreDict = df16.to_dict(orient='index')

#get standardized conservation scores
from sklearn import preprocessing
import numpy as np
import scipy.stats as stats
import itertools
with open('std_score.csv', 'w', encoding='UTF8') as f:
    writer=csv.writer(f)
    for key, value in NormScoreDict.items():
        keystr = [str(key)]
        row=[]
        vlist=[]
        for k, v in value.items():
            vlist.append(v)
        zscores = list(stats.zscore(vlist))
        clean_row = keystr+zscores
        #print(clean_row)
        writer.writerow(clean_row)
    print("done")

#Format normalized scores in dataframes
df16_2 = pd.read_csv('std_score.csv', delimiter=',', header=None)
df16_2.set_index([0], inplace=True)
df16_2.reset_index
df16_2.index.rename('Method', inplace=True)
df16_2.columns = pd.RangeIndex(0, len(df16.columns))

#Select variable residues
methods = ["KABAT,JORES,SCHNEIDER,SHENKIN,GERSTEIN,TAYLOR_GAPS,TAYLOR_NO_GAPS,ZVELIBIL,KARLIN,ARMON,THOMPSON,NOT_LANCET,MIRNY,WILLIAMSON,LANDGRAF,SANDER,VALDAR"]
print(methods)
cons_method = input("Please select which scoring method you would like to use: ")
score_strategy = input("Do you want to analyze using normalized or standardized scores? ")
if score_strategy == "normalized":
     scores = list(df16.loc["#{}".format(cons_method)])
     oc_set = set()
     pos = []
     res = []
     for idx, val in enumerate(scores):
         if val <= 0.3:
             if val not in oc_set:
                 oc_set.add(val)
                 pos.append(idx)
             else:
                 res.append(idx)
         else:
             continue
elif score_strategy == "standardized":
    scores = list(df16_2.loc["#{}".format(cons_method)])
    oc_set = set()
    pos = []
    res = []
    for idx, val in enumerate(scores):
        if val > 1:
            if val not in oc_set:
                oc_set.add(val)
                pos.append(idx)
            else:
                res.append(idx)
        else:
            continue
else:
    print("Error:Wrong input")

#Score dictionaries:
N_MethodDict = NormScoreDict["#{}".format(cons_method)]
StdScoreDict = df16_2.to_dict(orient='index')
S_MethodDict = StdScoreDict["#{}".format(cons_method)]

variable_pos = res + pos
variable_pos.sort()

var_list = []
for v in variable_pos:
    var_list.append(list(df12.loc[:, v]))
df22 = pd.DataFrame(var_list)
df23 = df22.T
df23.index = df12.index
actual_pos = []
for v in variable_pos:
    actual_pos.append(v + 1)
df23.columns = actual_pos

variable_res=list(df23.columns)

#Scrape consurf server
PATH = "/home/dell/WebDrivers/chromedriver"
chrome_options = webdriver.ChromeOptions()
preferences = {
    "profile.default_content_settings.popups": 0,
    "download.default_directory": os.getcwd() + os.path.sep,
    "directory_upgrade": True
}
chrome_options.add_experimental_option('prefs', preferences)
chrome_options.add_argument('--headless')
chrome_options.add_argument('window-size=1920x1080')
chrome_options.add_argument("disable-gpu")
s = Service(PATH)
driver = webdriver.Chrome(service=s, options=chrome_options)
driver.implicitly_wait(75)
driver.get("https://consurf.tau.ac.il/")
aa_button = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[2]/label/input")
aa_button.click()
no_structure = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[2]/div[2]/label/input")
no_structure.click()
yes_alignment = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[3]/div[1]/div[1]/label/input")
yes_alignment.click()
upload_aln = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[3]/div[2]/div[1]/input[1]")
alignment_file = os.path.abspath("sliced_seq_align.aln")
upload_aln.send_keys(alignment_file)
listofseq = []
for record in AlignIO.read("sliced_seq_align.aln", "clustal"):
    rec_id = (record.id)
    listofseq.append(rec_id)
print(listofseq)
seq_id = input("please type the sequence you want to analyze from the list given above: ")
try:
    ddmenu = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[3]/div[2]/div[1]/div[2]/select")
    ddmenu = Select(ddmenu)
    ddmenu.select_by_visible_text(str(seq_id))
except NoSuchElementException:
    seq_text = driver.find_element(by=By.NAME, value="msa_SEQNAME")
    seq_text.send_keys(str(listofseq[1]))
try:
    update_sel = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[3]/div[2]/div[1]/div[2]/a[2]")
    update_sel.click()
except NoSuchElementException:
    link_text = driver.find_element(by=By.LINK_TEXT, value="update selection")
    link_text.click()
try:
    no_tree = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[5]/div[1]/div[2]/label/input")
    no_tree.click()
except NoSuchElementException:
    tree_button = driver.find_element(by=By.NAME, value="TREE_yes_no")
    tree_button_no = tree_button.find_element(by=By.NAME, value="yes")
    tree_button_no.click()
try:
    no_email = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[6]/div[9]/label/input")
    no_email.click()
except NoSuchElementException:
    no_mail_check = driver.find_element(by=By.NAME, value="send_user_mail")
    no_mail_check.click()
try:
    submit_button = driver.find_element(by=By.XPATH, value="/html/body/div[3]/div/form/div[6]/input")
    submit_button.click()
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


def waitUntilDownloadCompleted(maxTime=600):
    driver.execute_script("window.open()")
    # switch to new tab
    driver.switch_to.window(driver.window_handles[-1])
    # navigate to chrome downloads
    driver.get('chrome://downloads')
    # define the endTime
    endTime = time.time() + maxTime
    while True:
        try:
            # get the download percentage
            downloadPercentage = driver.execute_script(
                "return document.querySelector('downloads-manager').shadowRoot.querySelector('#downloadsList downloads-item').shadowRoot.querySelector('#progress').value")
            # check if downloadPercentage is 100 (otherwise the script will keep waiting)
            if downloadPercentage == 100:
                # exit the method once it's completed
                return downloadPercentage
        except:
            pass
        # wait for 1 second before checking the percentage next time
        time.sleep(1)
        # exit method if the download not completed with in MaxTime.
        if time.time() > endTime:
            break
waitUntilDownloadCompleted(150)

print("downloaded consurf results.zip")

#Extract items from consurf_output.zip
import glob
from pathlib import Path
for name in glob.glob('Consurf*.zip'):
    zip_file = os.path.abspath(name)
    zip_path = Path(zip_file)
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
with zipfile.ZipFile(zip_path, "r") as zip_ref:
    zip_ref.extractall("/home/dell/Desktop/GTPase Algorithms/venv/consurf_output")

#Editing consurf grade files
lines = []
with open("/home/dell/Desktop/GTPase Algorithms/venv/consurf_output/consurf.grades", 'r') as fp:
    lines = fp.readlines()
with open("/home/dell/Desktop/GTPase Algorithms/venv/consurf_output/consurf_grades1", 'w') as fp:
    for number, line in enumerate(lines):
        if number > 13:
            fp.write(line)

#Read consurf grades using pandas dataframes
df24 = pd.DataFrame
with open('/home/dell/Desktop/GTPase Algorithms/venv/consurf_output/consurf_grades1', newline='') as grades:
    grade_reader = csv.reader(grades, delimiter='\t')
    list1 = []
    for grade in grade_reader:
        list1.append(grade)
df24 = list1
list2 = []
for data in list1:
    list3 = []
    for d in data:
        if d != '':
            list3.append(d)
    list2.append(list3)
df25 = pd.DataFrame(list2)
df25.fillna('', inplace=True)
df25.columns = [df25.iloc[0], df25.iloc[1]]
df25 = df25[2:-4]
df25.reset_index(drop=True, inplace=True)
df25 = df25.T.reset_index(level=1,drop=True).T
df25.columns = [i.strip() for i in df25.columns]
df25 = df25.applymap(lambda x : x.strip())
df25.POS = df25.POS.apply(int)
df26=df25[df25.POS.isin(list(df23.columns))][['POS','B/E','FUNCTION']].reset_index(drop=True)
exposed = list(df26[df26['B/E']=='e']['POS'])
buried = list(df26[df26['B/E']=='b']['POS'])
functional = list(df26[df26['FUNCTION']=='f']['POS'])
print("The exposed non-conserved residues are: {}".format(exposed))
print("The buried non-conserved residues are: {}".format(buried))
print("The functional non-conserved residues are: {}".format(functional))

#Scrape AlphaFold Database to get PDB structure of the required sequence
PATH = "/home/dell/WebDrivers/chromedriver"
chrome_options = webdriver.ChromeOptions()
preferences = {
                "profile.default_content_settings.popups": 0,
                "download.default_directory": '/home/dell/Desktop/GTPase Algorithms/venv/pdb_files/',
                "directory_upgrade": True
            }
chrome_options.add_experimental_option('prefs', preferences)
chrome_options.add_argument('--headless')
chrome_options.add_argument('window-size=1920x1080')
chrome_options.add_argument("disable-gpu")
s = Service(PATH)
driver=webdriver.Chrome(service=s, options=chrome_options)
driver.implicitly_wait(75)
driver.get("https://alphafold.ebi.ac.uk/")
listofseq=[]
for record in SeqIO.parse("sliced_sequences.fasta", "fasta"):
    listofseq.append(record.id)
print(listofseq)
seq_id=input("please select the sequence whose stucture you want to fetch from PDB: ")
seq_id_1 = seq_id[3:9]
print(seq_id_1)
driver.implicitly_wait(45)
query_input = driver.find_element(by=By.XPATH, value="/html/body/div[1]/app-root/app-header/div/div[2]/div/div/div[2]/app-header-search/div[1]/input").send_keys(seq_id_1)
search=driver.find_element(by=By.XPATH, value="/html/body/div[1]/app-root/app-header/div/div[2]/div/div/div[2]/app-header-search/button/span").click()
select_result= driver.find_element(by=By.CLASS_NAME, value="vf-link").click()
download_pdb=driver.find_element(by=By.LINK_TEXT, value="PDB file").click()
def waitUntilDownloadCompleted(maxTime=600):
    driver.execute_script("window.open()")
    # switch to new tab
    driver.switch_to.window(driver.window_handles[-1])
    # navigate to chrome downloads
    driver.get('chrome://downloads')
    # define the endTime
    endTime = time.time() + maxTime
    while True:
        try:
            # get the download percentage
            downloadPercentage = driver.execute_script(
                "return document.querySelector('downloads-manager').shadowRoot.querySelector('#downloadsList downloads-item').shadowRoot.querySelector('#progress').value")
            # check if downloadPercentage is 100 (otherwise the script will keep waiting)
            if downloadPercentage == 100:
                # exit the method once it's completed
                return downloadPercentage
        except:
            pass
        # wait for 1 second before checking the percentage next time
        time.sleep(1)
        # exit method if the download not completed with in MaxTime.
        if time.time() > endTime:
            break
waitUntilDownloadCompleted(75)
print("Downloaded PDB File")
def newest(path):
    files = os.listdir(path)
    paths = [os.path.join(path, basename) for basename in files]
    return max(paths, key=os.path.getctime)
PDB_Path = newest("/home/dell/Desktop/GTPase Algorithms/venv/pdb_files/")
os.rename(PDB_Path, r'/home/dell/Desktop/GTPase Algorithms/venv/pdb_files/{}_AFM.pdb'.format(seq_id_1))
PDB_Path1 = newest("/home/dell/Desktop/GTPase Algorithms/venv/pdb_files/")
PDB_name = os.path.basename(PDB_Path1)
print(PDB_name)

#Modify PDB File to add Conservation scores to B-Factor column
seq_index = seq_id[20:]
seq_index_list = seq_index.split("-")
seq_index_start = int(seq_index_list[0])
seq_index_stop = int(seq_index_list[1])
#change cwd
os.chdir('/home/dell/Desktop/GTPase Algorithms/venv/pdb_files/')
print(os.getcwd)
#slice pdb file into the respective selected domain
seleres = "python pdb_selres.py -{}:{} {} > sliced_{}".format(seq_index_start, seq_index_stop, PDB_name, PDB_name)
subprocess.call(seleres, shell=True)
print("done")

if score_strategy == "standardized":
    scores = S_MethodDict
elif score_strategy == "normalized":
    scores = N_MethodDict
else:
    print("wrong input")
values=[value for key, value in scores.items()]
keys = []
for i in range (seq_index_start,(seq_index_stop+1)):
    keys.append(i)
    #print (keys, values)
    range_scores = dict(zip(keys, values))
print(range_scores)

db1 = pdb2sql("sliced_{}".format(PDB_name))
#db1.print()
resseq = []
for x in db1.get('resSeq'):
    for i in x:
        resseq.append(i)
score_l = [range_scores[i] for i in resseq]
val = np.array(score_l, dtype=float)
db1.update_column('temp', values=val,)
#db1.print()
db1.exportpdb("temp_modified.pdb")
