#!/usr/bin/env python
# coding: utf-8

# In[ ]:

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

def Fetch_PFAM(fname):
    fin = []
    for record in SeqIO.parse(fname+".fasta", "fasta"):
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
        driver.find_element(by=By.XPATH, value='/html/body/div[5]/div[4]/div/a[1]').click()
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
    return(fin)

# In[ ]:


def slicedseq_PFAMdata(fin):
    df3 = pd.DataFrame()
    for f in fin:
        df = f[0]
        df2 = df.T.reset_index(level=1).T
        df2.columns = [i + "-" + j if (i != j and j != "") else i for i, j in zip(df2.columns, df2.loc['level_1'])]
        df2 = df2.drop('level_1', axis=0)
        df3 = pd.concat([df3, df2])
    data = []
    for record in SeqIO.parse(name + ".fasta", "fasta"):
        data.append([">" + record.id, "".join(record.seq)])
    df_seq = pd.DataFrame(data, columns=['id', 'seq'])
    df3 = df3.merge(df_seq, on='id')
    df3['sliced_sequence'] = df3.apply(lambda x: x['seq'][x['Envelope-Start'] - 1:x['Envelope-End']], axis=1)
    df3['new_id'] = df3.apply(lambda x: x['id'] + "|" + str(x['Envelope-Start']) + "-" + str(x['Envelope-End']), axis=1)
    return df3


# In[ ]:


def domain_dataframe(familyID, df3):
    ln = df3['Family'].count()
    df4 = pd.DataFrame()
    for x in range(ln):
        if df3['Family'].iloc[x].find(familyID) == 0:
            df4 = pd.concat([df4, df3.iloc[x:x + 1, :]])
        else:
            continue
    df4.reset_index(drop=True)
    return df4


# In[ ]:


def sliceddomains_FASTA(df4):
    with open("domain_seq.fasta", "w") as file:
        for idx, row in df4.iterrows():
            file.write(row['new_id'] + "\n")
            file.write(row['sliced_domainseqs'] + "\n")
    print("wrote sliced sequences in domain_seq.fasta")


# In[ ]:


def align_domains():
    align = "./probcons/probcons -clustalw domain_seq.fasta > domain_seq_aligned.aln"
    subprocess.call(align, shell=True)
    print ("Aligned domain sequences.")


# In[ ]:


def format_alignout():
    with open('domain_seq_aligned.aln', 'r') as file:
        data = file.readlines()
        data[0] = "CLUSTAL W (1.82) multiple sequence alignment\n"
    with open('sliced_seq_align.aln', 'w') as file:
        file.writelines(data)


# In[ ]:


def conservation_scores_calc():
    cons_score = "java -jar compbio-conservation-1.1.jar -i=domain_seq_aligned.aln -o=cons_out.csv"
    norm_score = "java -jar compbio-conservation-1.1.jar -i=domain_seq_aligned.aln -n -o=norm_out.csv"
    subprocess.call(norm_score, shell=True)
    subprocess.call(cons_score, shell=True)
    print("calculated conservation scores for the alignment")


# In[ ]:


def conservation_scores_df():
    df5 = pd.read_csv('cons_out.csv', delimiter='\s+', header=None)
    df5.set_index([0], inplace=True)
    df5.reset_index
    df5.index.rename('Method', inplace=True)
    df5.columns = pd.RangeIndex(0, len(df5.columns))
    df5.to_csv('conservation_scores.csv')
    return df5


# In[ ]:


def normalized_scores_df():
    df6 = pd.read_csv('norm_out.csv', delimiter='\s+', header=None)
    df6.set_index([0], inplace=True)
    df6.reset_index
    df6.index.rename('Method', inplace=True)
    df6.columns = pd.RangeIndex(0, len(df16.columns))
    df6.to_csv('normalized_scores.csv')
    NormScoreDict = df6.to_dict(orient='index')
    return df6, NormScoreDict


# In[ ]:


def predicted_variable_res(df6):
    methods = ["KABAT,JORES,SCHNEIDER,SHENKIN,GERSTEIN,TAYLOR_GAPS,TAYLOR_NO_GAPS,ZVELIBIL,KARLIN,ARMON,THOMPSON,NOT_LANCET,MIRNY,WILLIAMSON,LANDGRAF,SANDER,VALDAR"]
    print(methods)
    cons_method = input("Please select which scoring method you would like to use: ")
    scores = list(df6.loc["#{}".format(cons_method)])
    oc_set = set()
    pos = []
    res = []
    for idx, val in enumerate(scores):
        if val <= 0.2:
            if val not in oc_set:
                oc_set.add(val)
                pos.append(idx)
            else:
                res.append(idx)
    variable_pos = res + pos
    variable_pos.sort()
    print("The predicted variable and non-conserved residues are: ", variable_pos)
    return variable_pos


# In[ ]:


def variable_res_df():
    alignment = AlignIO.read("domain_seq_aligned.aln", "clustal")
    l = []
    for i in alignment.__dict__['_records']:
        l.append([i.id, "".join(i.seq)])
    df7 = pd.DataFrame(l)
    df7.columns = ['id', 'seq']
    df8 = pd.DataFrame([list(word) for word in list(df7.seq)])
    df9 = pd.DataFrame([list(word) for word in list(df4.sliced_sequence)], index=df4.index)
    var_list = []
    for v in variable_pos:
        var_list.append(list(df9.loc[:, v]))
    df10 = pd.DataFrame(var_list)
    df11 = df10.T
    df11.index = df4.index
    actual_pos = []
    for v in variable_pos:
        actual_pos.append(v)
    df11.columns = actual_pos
    return df11


# In[ ]:


def protein_domain_varible_pos():
    df12 = pd.DataFrame([list(word) for word in list(df7.seq)], index=df7.index)
    df13 = []
    for idx,row in df11.iterrows():
        if idx!="Consensus":
            offset = int(idx.split("|")[-1].split("-")[0])
            s = list(df12.loc[idx,:])
            s1 = list(enumerate(s))
            s2 = [i for i in s1 if i[1]!="-"]
            s3 = [(i[0],i[1],idx+offset) for idx,i in enumerate(s2)]
            idxs = list(df23.columns)
            df13.append([j[2] for i in idxs for j in s3 if j[0]==i])
        else:
            df13.append([])
    df14 = pd.DataFrame(df13,columns= list(df23.columns))
    df14['id'] = list(df11.index)
    df14 = df14.set_index('id')
    return df14


# In[ ]:


def Consurf_Scraper():
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
    alignment_file = os.path.abspath("domain_seq_aligned.aln")
    upload_aln.send_keys(alignment_file)
    listofseq = []
    for record in AlignIO.read("domain_seq_aligned.aln", "clustal"):
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


# In[ ]:


def waitUntilDownloadCompleted(maxTime=600):
    driver.execute_script("window.open()")
    driver.switch_to.window(driver.window_handles[-1])
    driver.get('chrome://downloads')
    endTime = time.time() + maxTime
    while True:
        try:
            downloadPercentage = driver.execute_script(
                "return document.querySelector('downloads-manager').shadowRoot.querySelector('#downloadsList downloads-item').shadowRoot.querySelector('#progress').value")
            if downloadPercentage == 100:
                return downloadPercentage
        except:
            pass
        time.sleep(1)
        if time.time() > endTime:
            break


# In[ ]:


def extract_consurf_results():
    import glob
    from pathlib import Path
    for name in glob.glob('Consurf*.zip'):
        zip_file = os.path.abspath(name)
        zip_path = Path(zip_file)
        cwd = os.getcwd()
    new_folder_name = outputs
    output_path = os.path.join(cwd,new_folder_name)
    os.mkdir(output_path)
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(output_path)
    return output_path


# In[ ]:


def modify_consurfgrades ():
    lines = []
    with open(output_path+"/consurf_grades", 'r') as fp:
        lines = fp.readlines()
    with open(output_path+"/consurf_grades1", 'w') as fp:
        for number, line in enumerate(lines):
            if number > 13:
                fp.write(line)


# In[ ]:


def consurfgrades_df ():
    df15 = pd.DataFrame
    with open(output_path+"/consurf_grades1", newline='') as grades:
        grade_reader = csv.reader(grades, delimiter='\t')
        list1 = []
        for grade in grade_reader:
            list1.append(grade)
    df15 = list1
    list2 = []
    for data in list1:
        list3 = []
        for d in data:
            if d != '':
                list3.append(d)
        list2.append(list3)
    df16 = pd.DataFrame(list2)
    df16.fillna('', inplace=True)
    df16.columns = [df16.iloc[0], df16.iloc[1]]
    df16 = df16[2:-4]
    df16.reset_index(drop=True, inplace=True)
    df16 = df16.T.reset_index(level=1,drop=True).T
    df16.columns = [i.strip() for i in df16.columns]
    df16 = df16.applymap(lambda x : x.strip())
    df16.POS = df16.POS.apply(int)
    df17=df16[df16.POS.isin(list(df14.loc[seq_id,:]))][['POS','B/E','FUNCTION']].reset_index(drop=True)
    exposed = list(df17[df17['B/E']=='e']['POS'])
    buried = list(df17[df17['B/E']=='b']['POS'])
    functional = list(df17[df17['FUNCTION']=='f']['POS'])
    print("The exposed non-conserved residues are: {}".format(exposed))
    print("The buried non-conserved residues are: {}".format(buried))
    print("The functional non-conserved residues are: {}".format(functional))
    df17.to_csv("variable_res.csv")
    return df17


# In[ ]:


def create_results ():
    with open ('CPDA_Result.txt', 'w') as f:
        f.writelines('The variable residues based on the scoring strategy are: ')
        f.writelines( list( "%s\n" % item for item in variable_res) )
        f.writelines('The exposed variable residues are: ')
        f.writelines( list( "%s\n" % item for item in exposed))
        f.writelines('The buried variable residues are: ')
        f.writelines( list( "%s\n" % item for item in buried) )
        f.writelines('The functional variable residues are: ')
        f.writelines( list( "%s\n" % item for item in functional) )


# In[ ]:


def get_Alphafold_PDB():
    PATH = "/home/dell/WebDrivers/chromedriver"
    chrome_options = webdriver.ChromeOptions()
    preferences = {
                    "profile.default_content_settings.popups": 0,
                    "download.default_directory": output_path,
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
    return seq_id, seq_id_1


# In[ ]:


def waitUntilDownloadCompleted(maxTime=600):
    driver.execute_script("window.open()")
    driver.switch_to.window(driver.window_handles[-1])
    driver.get('chrome://downloads')
    endTime = time.time() + maxTime
    while True:
        try:
            downloadPercentage = driver.execute_script(
                "return document.querySelector('downloads-manager').shadowRoot.querySelector('#downloadsList downloads-item').shadowRoot.querySelector('#progress').value")
            if downloadPercentage == 100:
                return downloadPercentage
        except:
            pass
        time.sleep(1)
        if time.time() > endTime:
            print("file not downloaded")
            break


# In[ ]:


def newest(path):
    files = os.listdir(path)
    paths = [os.path.join(path, basename) for basename in files]
    return max(paths, key=os.path.getctime)


# In[ ]:


def PDB_modifier ():
    seq_index = seq_id[20:]
    seq_index_list = seq_index.split("-")
    seq_index_start = int(seq_index_list[0])
    seq_index_stop = int(seq_index_list[1]) 
    os.chdir(output_path)
    seleres = "python pdb_selres.py -{}:{} {} > sliced_{}".format(seq_index_start, seq_index_stop, PDB_name, PDB_name)
    subprocess.call(seleres, shell=True)
    print("done")


# In[ ]:


def Bfactore_modifier ():
    scores = N_MethodDict
    print("wrong input")
    values=[value for key, value in scores.items()]
    keys = []
    for i in range (seq_index_start,(seq_index_stop+1)):
        keys.append(i)
        #print (keys, values)
        range_scores = dict(zip(keys, values))
    modified_pdb= ("sliced_{}".format(PDB_name))
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
    return modified_pdb, range_scores

