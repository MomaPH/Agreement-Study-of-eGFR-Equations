#CKD-EPI 2021 function source https://www.kidney.org/ckd-epi-creatinine-equation-2021-0
def CKD_EPI(Sex, Age, Scr):
  if Sex == 'M':
    k = 0.9
    a=-0.302
    r= Scr/k
    if r<=1:
      eGFR=142*(r**a)*(0.9938**Age)
    elif r>1:
      eGFR=142*(r**(-1.200))*(0.9938**Age)
  elif Sex == 'F':
    k = 0.7
    a=-0.241
    r= Scr/k
    if r<=1:
      eGFR=142*(r**a)*(0.9938**Age)*1.012
    elif r>1:
      eGFR=142*(r**(-1.200))*(0.9938**Age)*1.012
  return round(eGFR, 1)

#MDRD without race factor function source https://www.niddk.nih.gov/research-funding/research-programs/kidney-clinical-research-epidemiology/laboratory/glomerular-filtration-rate-equations/adults/previous#mdrd
def I4MDRD(Sex, Age, Scr):
 MeGFR=175*(Scr**-1.154)*(Age**-0.203)
 if Sex=='F':
    MeGFR*=0.742
 return round(MeGFR, 1)

#EKFC function source DOI:10.7326/M20-4366
def EKFC(Sex, Age, Scr):
  #Importing the math module to compute ln and exponential functions
  import math
  #calculating Q for different scenarios of age and sex
  if Age>25 and Sex=='M':
    Q=0.90
    r=Scr/Q
    #further tuning the formula for different age ranges
    if Age<=40:
      #finetuning for different r
      if r<1:
        EeGFr=107.3*(r**(-0.322))
      elif r>=1:
        EeGFr=107.3*(r**-1.132)
    elif Age>40:
      if r<1:
        EeGFr=107.3*(r**-0.322)*0.990**(Age-40)
      elif r>=1:
        EeGFr=107.3*(r**-1.132)*0.990**(Age-40)
  elif Age>25 and Sex=='F':
    Q=0.70
    r=Scr/Q
    if Age<=40:
      #finetuning for different r
      if r<1:
        EeGFr=107.3*(r**(-0.322))
      elif r>=1:
        EeGFr=107.3*(r**-1.132)
    elif Age>40:
      if r<1:
        EeGFr=107.3*(r**-0.322)*0.990**(Age-40)
      elif r>=1:
        EeGFr=107.3*(r**-1.132)*0.990**(Age-40)
  elif Age<=25 and Sex=='M':
    Q=math.exp( 3.200 + 0.259*Age-0.543*math.log(Age) - 0.00763*Age**2 + 0.0000790*Age**3)
    r=Scr/Q
    if r<1:
      EeGFr=107.3*(r**-0.322)
    elif r>=1:
      EeGFr=107.3*(r**-1.132)
  elif Age<=25 and Sex=='F':
    Q=math.exp(3.080 + 0.177*Age-0.223*math.log(Age) - 0.00596*Age**2 + 0.0000686*Age**3)
    r=Scr/Q
    if r<1:
      EeGFr=107.3*(r**-0.322)
    elif r>=1:
      EeGFr=107.3*(r**-1.132)
  return round(EeGFr, 1)

# calculating random values to check with online calculators
print(f"MDRD= {I4MDRD('M', 60, 1.5)}, CKD-epi= {CKD_EPI('M', 60, 1.5)}, EKFC= {EKFC('M', 60, 1.5)}")
#check success

import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pyCompare as pyc
import seaborn as sbn
import pandas as pd
import scipy.stats as stats
import pingouin as pg
import numpy as np

vars1=pd.read_sas('/content/drive/MyDrive/BIOPRO_J.xpt')
vars1=pd.DataFrame(vars1)
col1keep=['SEQN','LBXSCR']
vars1=vars1[col1keep]
vars1.head(2)

vars2=pd.read_sas('/content/drive/MyDrive/DEMO_J20172018.xpt')
vars2=pd.DataFrame(vars2)
vars2.head(2)
col2keep=['SEQN','RIAGENDR','RIDAGEYR']
vars2=vars2[col2keep]
vars2.head(2)

vars=pd.merge(vars2, vars1, on='SEQN')
vars['RIAGENDR']= vars['RIAGENDR'].replace({1:'M', 2:'F'})
vars.dropna(inplace=True)
vars.head(2)

vars['egfr_mdrd'] = vars.apply(lambda row: I4MDRD(row['RIAGENDR'], row['RIDAGEYR'], row['LBXSCR']), axis=1)
vars['egfr_CKDEPI21'] = vars.apply(lambda row: CKD_EPI(row['RIAGENDR'], row['RIDAGEYR'], row['LBXSCR']), axis=1)
vars['egfr_EKFC'] = vars.apply(lambda row: EKFC(row['RIAGENDR'], row['RIDAGEYR'], row['LBXSCR']), axis=1)
vars.head()

vars = vars[vars["RIDAGEYR"] >= 18]
vars = vars[vars["RIDAGEYR"] <71]
vars=vars.head(200) #taking only 200 data points to avoid overloading the graphs
Desc=vars.describe()
Desc

# Create age bins for better visualization with violin plots
bins = [18, 24, 39, 64, 150]# Define age ranges
labels =["18-24", "25-39", "40-64", "65+"]
len(labels)==len(bins)-1 # Create labels for the bins
vars['age_bins'] = pd.cut(vars['RIDAGEYR'], bins=bins, labels=labels, right=False)
sbn.violinplot(data=vars, y='age_bins', x='egfr_mdrd', hue='age_bins', cut=0, density_norm='count')
plt.title('eGFR (MDRD) Distribution by Age Groups')
plt.ylabel('Age Group (years)')
plt.xlabel('eGFR (MDRD)')
#plt.tight_layout() # Adjust layout to prevent labels from being cut off
plt.show()

sbn.violinplot(data=vars, y='age_bins', x='egfr_CKDEPI21', hue='age_bins', density_norm='count', cut=0)
plt.title('eGFR (CKD-EPI) Distribution by Age Groups')
plt.ylabel("Age group (years)")
plt.xlabel("eGFR (CKD-EPI)")
plt.show()

sbn.violinplot(data=vars, x='egfr_EKFC', y='age_bins', hue='age_bins', density_norm='count', cut=0)
plt.title('eGFR (EKFC) Distribution by Age Groups')
plt.ylabel("Age group (years)")
plt.xlabel("eGFR (EKFC)")
plt.show()

fig, ax = plt.subplots(figsize=(8,5)) #choose fibonnaci suite number to obtain golden ratio and flex!!!
pyc.blandAltman(
    vars['egfr_mdrd'],
    vars['egfr_CKDEPI21'],
    title='MDRD & CKDEPI21 Simple Bland–Altman',
    loaColour='red',
    confidenceIntervalMethod='approximate',
    confidenceInterval=95,
    ax=ax
    )
ax.set_ylabel('DIFFERENCE BETWEEN MDRD & CKDEPI21')
ax.set_xlabel('MEAN OF MDRD & CKDEPI21')
plt.show()

plt.close('all') #avoid conflict with previous plot
fig, ax = plt.subplots(figsize=(8,5))
pyc.blandAltman(
    vars['egfr_mdrd'],
    vars['egfr_EKFC'],
    title='MDRD & EKFC Simple Bland–Altman',
    loaColour='red',
    confidenceIntervalMethod='approximate',
    confidenceInterval=95,
    ax=ax
    )
ax.set_ylabel('DIFFERENCE BETWEEN MDRD & EKFC')
ax.set_xlabel('MEAN OF MDRD & EKFC')
plt.show()

plt.close('all')
fig, ax = plt.subplots(figsize=(8,5))
pyc.blandAltman(
    vars['egfr_CKDEPI21'],
    vars['egfr_EKFC'],
    title='CKDEPI21 & EKFC Simple Bland–Altman',
    loaColour='red',
    confidenceIntervalMethod='approximate',
    confidenceInterval=95,
    ax=ax)
ax.set_ylabel('DIFFERENCE BETWEEN CKDEPI21 & EKFCRD')
ax.set_xlabel('MEAN OF KDEPI21 & EKFC')
plt.show()

#bland altman with percent difference
plt.close('all')
fig, ax = plt.subplots(figsize=(8,5))
pyc.blandAltman(
    vars['egfr_mdrd'],
    vars['egfr_CKDEPI21'],
    title='MDRD & CKDEPI21 Percentage Difference Bland–Altman',
    loaColour='red',
    percentage=True,
    confidenceIntervalMethod='approximate',
    confidenceInterval=95,
    ax=ax
    )
ax.set_ylabel('DIFFERENCE BETWEEN MDRD & CKDEPI21')
ax.set_xlabel('MEAN OF MDRD & CKDEPI21')
plt.show()

#bland altman with percent difference
plt.close('all')
fig, ax = plt.subplots(figsize=(8,5))
pyc.blandAltman(
    vars['egfr_mdrd'],
    vars['egfr_EKFC'],
    title='MDRD & EKFC Percentage Difference Bland–Altman',
    loaColour='red',
    confidenceIntervalMethod='approximate',
    percentage=True,
    confidenceInterval=95,
    ax=ax
    )
ax.set_ylabel('DIFFERENCE BETWEEN MDRD & EKFC')
ax.set_xlabel('MEAN OF MDRD & EKFC')
plt.show()

#bland altman with percent difference
plt.close('all')
fig, ax = plt.subplots(figsize=(8,5))
pyc.blandAltman(
    vars['egfr_CKDEPI21'],
    vars['egfr_EKFC'],
    title='CKDEPI21 & EKFC Percentage Difference Bland–Altman',
    loaColour='red',
    confidenceIntervalMethod='approximate',
    percentage=True,
    confidenceInterval=95,
    ax=ax)
ax.set_ylabel('DIFFERENCE BETWEEN CKDEPI21 & EKFC')
ax.set_xlabel('MEAN OF CKDEPI21 & EKFC')
plt.show()

Lowegfr=vars[vars["egfr_CKDEPI21"] <=90] #because all egfr equations in this study underestimate gfr compared to CKDEPI, we take the values where CKDEPI is inferior to the 90
Lowegfr.shape #checking if the number of data points is close to the recommendation of Bland altman

pyc.blandAltman(
    Lowegfr['egfr_mdrd'],
    Lowegfr['egfr_CKDEPI21'],
    title='LOW MDRD & CKDEPI21 Percentage Difference Bland–Altman',
    loaColour='red',
    confidenceIntervalMethod='approximate',
    confidenceInterval=95,
    percentage= True,
    )

plt.show()

pyc.blandAltman(
    Lowegfr['egfr_mdrd'],
    Lowegfr['egfr_EKFC'],
    title='low MDRD & EKFC percentage difference Bland–Altman',
    loaColour='red',
    confidenceIntervalMethod='approximate',
    confidenceInterval=95,
    percentage=True
        )
plt.show()

pyc.blandAltman(
    Lowegfr['egfr_CKDEPI21'],
    Lowegfr['egfr_EKFC'],
    title='low CKDEPI21 & EKFC Percentage Difference  Bland–Altman',
    loaColour='red',
    confidenceIntervalMethod='approximate',
    confidenceInterval=95,
    percentage=True
    )
plt.show()

