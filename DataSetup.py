""" Some script to prepare the raw data for analysis """

import pandas as pd
import numpy as np

# ======= Import the data =======

data = pd.read_csv("C:/Users/jamie/OneDrive/Documents/MMath Project/Summer2021/RFE002_demo_baseline.RCS.csv")
sampleSize = data.shape[0]

# ======= Some data formatting =======

# convert entries to dates
data['DateOfFirstConsultationMAT'] = pd.to_datetime(data['DateOfFirstConsultationMAT'], dayfirst=True)

# convert entries to numeric
data['outcome'] = pd.to_numeric(data['outcome'], errors='coerce')

# convert Yes/No data to 1-0 format
def yesnotoindicator(string):
    """ takes a {"Yes", "No"} valued variable and converts it to an integer.
        1 for "Yes", and 0 for "No" """
        
    if string=="Yes" or string=="No":
        return (string=="Yes")*1
    
    else:
        return string
    
def yesnoextoindicator(string):
    if string=="Yes":
        return 1
    elif string=="No":
        return 0
    else:
        return string

data['SmokingMAT'] = pd.to_numeric(data['SmokingMAT'].apply(yesnotoindicator, 1), errors='coerce')
data['IsSmoker'] = pd.to_numeric(data['SmokingMAT'].apply(yesnoextoindicator, 1), errors='coerce')
data['PolycysticOvaries'] = pd.to_numeric(data['PolycysticOvaries'].apply(yesnotoindicator, 1), errors='coerce')

# ======= Compute number of previous miscarriages =======

def naToZero(x):
    if pd.isna(x):
        return 0
    else:
        return x
    
data['PreSecondTrimesterMis'] = pd.to_numeric(data['PreSecondTrimesterMis'].apply(naToZero, 1), errors='coerce')
data['PreClinicMiscarriages'] = data['PreFirstTrimesterMis'] + data['PreSecondTrimesterMis']

# ======= Create Success Indicator for First Pregnancy =======

def isPregSuccess(x):
    if pd.isna(x) or x==10:
        return pd.NA
    elif x==1:
        return 1
    else:
        return 0

data['FirstPregSuccess'] = data['outcome'].apply(isPregSuccess, 1)
    
# ======= Fill in Missing Coffee Values =======

def emptyToZero(x):
    if pd.isna(x):
        return 0
    else:
        return x
    
data['Coffee'] = data['CoffeeWeekMAT'].apply(emptyToZero, 1)
data['Tea'] = data['TeaWeekMAT'].apply(emptyToZero, 1)
data['Caffeine'] = data['Coffee'] + 0.5*data['Tea']

# ======= Fill in Missing Time to Event Data =======

# A substantial amount of women have missing data for the time to conception/ no event
# for many women this value is recoverable from other data
# =============================================================================
# for i in range(sampleSize):
#     if np.isnan(data['Time to conception/no event'][i]):
#         # if no conception is recorded and we can, borrow the value from Time to viable pregnancy/no event 
#         if data['Conception'][i]==0 & np.logical_not(np.isnan(data['Time to viable pregnancy/no event'][i])):
#             data['Time to conception/no event'][i] = data['Time to viable pregnancy/no event'][i]
#             
#         # if we can't, unfortunately the date of extraction is not consistent, so we are stuck :(
#         # Thankfully, a manual inspection reveals that this only affects 3 data points :)
#         
#         # If a conception is recorded and it is successful, we can subtract 24 weeks (168 days) off
#         # of Time to viable pregnancy/no event        
#         elif (data['Conception'][i]==1) & (data['outcome'][i]==1):
#             if np.logical_not(np.isnan(data['Time to viable pregnancy/no event'][i])):
#                 data['Time to conception/no event'][i] = np.maximum(0, data['Time to viable pregnancy/no event'][i] - 168)
# =============================================================================
  
extractiondate = pd.to_datetime('01/07/2021', dayfirst=True)
v = extractiondate - data['DateOfFirstConsultationMAT']
data['dayssincefirstconsultation'] = v.dt.days

for i in range(sampleSize):
    if np.isnan(data['Time to conception/no event'][i]):
        # if no conception is recorded and if we can, we can borrow the value from Time to viable pregnancy/no event
        # OR: compute the noumber of days between first consultation and data extraction
        if data['Conception'][i]==0:
            data['Time to conception/no event'][i] = data['dayssincefirstconsultation'][i]
            
        # if we can't, unfortunately the date of extraction is not consistent, so we are stuck :(
        # Thankfully, a manual inspection reveals that this only affects 3 data points :)
        
        # If a conception is recorded and it is successful, we can subtract 24 weeks (168 days) off
        # of Time to viable pregnancy/no event        
        elif (data['Conception'][i]==1) & (data['outcome'][i]==1):
            if np.logical_not(np.isnan(data['Time to viable pregnancy/no event'][i])):
                data['Time to conception/no event'][i] = np.maximum(0, data['Time to viable pregnancy/no event'][i] - 168)
    
        
        




