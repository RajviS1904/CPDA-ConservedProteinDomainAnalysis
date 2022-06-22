import requests
import os
from bs4 import BeautifulSoup
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
import csv

pdb_paths=['/home/dell/Desktop/CPDAS/venv/pdb-files/P0A6P5_4_119.pdb', '/home/dell/Desktop/CPDAS/venv/pdb-files/P0A6P5_204_322.pdb', '/home/dell/Desktop/CPDAS/venv/pdb-files/P50743_5_120.pdb', '/home/dell/Desktop/CPDAS/venv/pdb-files/P50743_177_295.pdb']
PATH = "/home/dell/WebDrivers/chromedriver"
chrome_options = webdriver.ChromeOptions()
preferences = {
    "profile.default_content_settings.popups": 0,
    "download.default_directory": os.getcwd() + os.path.sep,
    "directory_upgrade": True
}
s = Service(PATH)
driver = webdriver.Chrome(service=s, options=chrome_options)
for path in pdb_paths:
    email = "rajvisrivastav@ggmail.com" #input('Please enter your email ID')
    driver.implicitly_wait(20)
    driver.get("http://curie.utmb.edu/getarea.html")
    pdb_upload = driver.find_element(by=By.XPATH, value="/html/body/center/form/table/tbody/tr[1]/td/div/large/input[1]").send_keys(path)
    email = driver.find_element(by=By.XPATH, value='/html/body/center/form/table/tbody/tr[1]/td/div/large/input[4]').send_keys(email)
    submit = driver.find_element(by=By.XPATH, value='/html/body/center/form/table/tbody/tr[1]/td/div/large/input[5]').click()
    WebDriverWait(driver, 175).until(EC.presence_of_element_located((By.CSS_SELECTOR, "head > title")))
    soup = BeautifulSoup(driver.page_source, "lxml")
    prelist = soup.find_all('pre')
    info = []
    for pre in prelist:
        info.append(pre.text)
    file_name = (os.path.basename(path)).split('.')
    with open ('getarea.txt', 'w') as file:
            for i in info:
                file.write('\n'+ i)
    with open ('getarea.txt') as fin, open( file_name[0]+'.csv', 'w') as fout:
        o=csv.writer(fout)
        for line in fin:
            o.writerow(line.split())
print('done')
