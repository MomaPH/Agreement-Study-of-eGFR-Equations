# Agreement Study of eGFR Equations (NHANES 2017–2018)

Compare agreement between **EKFC**, **CKD-EPI 2021 (race-free, IDMS-traceable)**, and **4V-MDRD (IDMS)** using **NHANES 2017–2018** data. Emphasis is on **clinical decision thresholds**, specifically the impact on:

- **CKD staging reclassification**
- **drug dosing thresholds** (initial focus: **metformin** and **tenofovir**)

> This repository is intended as a quantitative, reproducible agreement study with clinically interpretable outputs.

---

## Study design

### Data source
- NHANES 2017–2018 (CDC)
- Key variables (current):
  - Serum creatinine: **LBXSCR** (mg/dL)
  - Age: **RIDAGEYR** (years)
  - Sex: **RIAGENDR**
  
### Inclusion criteria (current)
- Adults **18–70 years** (RIDAGEYR >= 18 and < 71)
- First **200 rows** used temporarily for plotting/performance

---

## eGFR equations included

1. **CKD-EPI 2021** (race-free, creatinine-based, IDMS)
2. **4-variable MDRD** (IDMS-traceable form)
3. **EKFC**

Implementation details and references are documented in `/src/egfr_equations.py` 

---

## Primary outcomes

### Agreement
Implemented:
- Bland–Altman analysis (bias + limits of agreement)
  
Planned additions:
- ICC
- Alluvial plots for CKD reclassifciation and 
- CKD stage reclassification table
- Drug dosing reclassification at thresholds relevant to:
  - **Metformin**
  - **Tenofovir**
- Cohen’s kappa
- Gwet’s AC1 (recommended when prevalence/imbalance may distort kappa)

---

## Repository structure

