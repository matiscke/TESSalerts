#This script add the columns (Teff, magnitudes,log_g, etc) from the TICs catalog-tables to the Alerts_TESS_Table (targets selected only if are visible from Calar Alto)

from astropy.table import Table
import numpy as np
import os
import pandas as pd
import gc

gc.enable()


def search_closest_even(n):
    int_n = int(n)
    if int_n%2 == 0:
        if int_n>0:
            return (int_n,int_n+2)
        elif int_n<0:
            return (int_n-2,int_n)
        else:
            if n>=0:
                return (int_n,int_n+2)
            else:
                return (int_n-2,int_n)
    else:
        return (int_n-1,int_n+1)


names = np.linspace(1,89,89).astype(int).astype(str) 
n = len(names)
for i in range(n):
    names[i] = 'col'+names[i]

#The table Alerts_TESS_Table.dat are the table-output from tess_alerts.py
t_alerts = Table.read('Alerts_TESS_Table.dat',format = 'ascii')

tic_id = t_alerts['tic_id'].data.data
dec_alerts = t_alerts['DEC'].data.data


n = len(t_alerts)
col_table_tic_name = []
for i in range(n):
    #print(str(i)+'/'+str(n))
    if dec_alerts[i] !='...' and tic_id[i]!='...':
        dec = float(dec_alerts[i])
        limits_dec = np.array(search_closest_even(dec))
        letras = np.where(limits_dec>=0,'N','S')
        limits_dec = np.abs(limits_dec).astype(str)
        limits_dec[0] = limits_dec[0].zfill(2)
        limits_dec[1] = limits_dec[1].zfill(2)
        table_tic_name = 'tic_dec'+limits_dec[0] +'_00'+letras[0] +'__'+ limits_dec[1]+'_00'+letras[1] +'.csv'
        col_table_tic_name.append(table_tic_name)
        

col_table_tic_name = np.array(col_table_tic_name)      
table_tic_uniq = np.unique(col_table_tic_name)
table_tic_uniq = table_tic_uniq[np.where(table_tic_uniq!='...')]

col_Teff = ['...']*len(t_alerts)
col_DTeff = ['...']*len(t_alerts)
col_gaia = ['...']*len(t_alerts)
col_Jmag = ['...']*len(t_alerts)
col_eJmag = ['...']*len(t_alerts)
col_Hmag = ['...']*len(t_alerts)
col_eHmag = ['...']*len(t_alerts)
col_Kmag = ['...']*len(t_alerts)
col_eKmag = ['...']*len(t_alerts)
col_logg = ['...']*len(t_alerts)
col_elogg = ['...']*len(t_alerts)



m = len(table_tic_uniq)
count = 0
for element in table_tic_uniq:
    print(str(count)+'/'+str(m))
    count+=1
    k = np.where(col_table_tic_name == element)[0]
    if len(k)>0:
        #t_tess = pd.read_csv(element,names = names)
        for t_tess in pd.read_csv(element,chunksize=1E6,names=names):
            tic_id_tess = t_tess['col1'].get_values()
            tess_teff = t_tess['col65'].get_values()
            tess_dteff = t_tess['col66'].get_values()
            tess_logg = t_tess['col67'].get_values()
            tess_elogg = t_tess['col68'].get_values()
            tess_gaia = t_tess['col9'].get_values()           
            tess_Jmag = t_tess['col43'].get_values()
            tess_eJmag = t_tess['col44'].get_values()
            tess_Hmag = t_tess['col45'].get_values()
            tess_eHmag = t_tess['col46'].get_values()
            tess_Kmag = t_tess['col47'].get_values()
            tess_eKmag = t_tess['col48'].get_values()
            for i in k:
                t_alerts_tic = int(tic_id[i])
                j = np.where(tic_id_tess==t_alerts_tic)[0]
                if len(j)>0:
                    j = j[0]
                    teff = tess_teff[j]
                    dteff = tess_dteff[j]
                    gaia = tess_gaia[j]
                    logg = tess_logg[j]
                    elogg = tess_elogg[j]                  
                    Jmag = tess_Jmag[j]
                    eJmag = tess_eJmag[j]
                    Hmag = tess_Hmag[j]
                    eHmag = tess_eHmag[j]
                    Kmag = tess_Kmag[j]
                    eKmag = tess_eKmag[j]
                    if ~np.isnan(teff):
                        col_Teff[i] = teff
                    if ~np.isnan(dteff):
                        col_DTeff[i] = dteff
                    if ~np.isnan(logg):
                        col_logg[i] = logg
                    if ~np.isnan(elogg):
                        col_elogg[i] = elogg                        
                    if ~np.isnan(gaia):
                        col_gaia[i] = gaia
                    if ~np.isnan(Jmag):
                        col_Jmag[i] = Jmag
                    if ~np.isnan(eJmag):
                        col_eJmag[i] = eJmag                        
                    if ~np.isnan(Hmag):
                        col_Hmag[i] = Hmag
                    if ~np.isnan(eHmag):
                        col_eHmag[i] = eHmag                       
                    if ~np.isnan(Kmag):
                        col_Kmag[i] = Kmag
                    if ~np.isnan(eKmag):
                        col_eKmag[i] = eKmag                        

     
from astropy.table import Column
col1 = Column(data = col_gaia,name = 'GAIA_ID', dtype = str)
col2 = Column(data = col_Teff,name = 'Teff', dtype = str)
col3 = Column(data = col_DTeff,name = 'e_Teff', dtype = str)
col4 = Column(data = col_logg,name = 'log_g', dtype = str)
col5 = Column(data = col_elogg,name = 'e_log_g', dtype = str)
col6 = Column(data = col_Jmag,name = 'Jmag', dtype = str)
col7 = Column(data = col_eJmag,name = 'e_Jmag', dtype = str)
col8 = Column(data = col_Hmag,name = 'Hmag', dtype = str)
col9 = Column(data = col_eHmag,name = 'e_Hmag', dtype = str)
col10 = Column(data = col_Kmag,name = 'Kmag', dtype = str)
col11= Column(data = col_eKmag,name = 'e_Kmag', dtype = str)

t_alerts.add_column(col1)
t_alerts.add_column(col2)
t_alerts.add_column(col3)
t_alerts.add_column(col4)
t_alerts.add_column(col5)
t_alerts.add_column(col6)
t_alerts.add_column(col7)
t_alerts.add_column(col8)
t_alerts.add_column(col9)
t_alerts.add_column(col10)
t_alerts.add_column(col11)

#Table with all targets visible from Calar Alto with columns (Gaia_ID, Teff, log_g, etc...)
t_alerts.write('Alerts_TESS_Table_columns.csv',format = 'csv',overwrite=True)
t_alerts.write('Alerts_TESS_Table_columns.dat',format = 'ascii',overwrite=True)

#select only the TIC targets with Teff<4000K
t_alerts = Table.read('Alerts_TESS_Table_columns.dat',format = 'ascii')
tic_teff = t_alerts['Teff'].data
final_teff=np.where(tic_teff<4000)[0] 
selected=t_alerts[final_teff] 

#Final Alerts table
selected.write('Final_Alerts_TESS_Table_columns.csv',format = 'csv',overwrite=True)
selected.write('Final_Alerts_TESS_Table_columns.dat',format = 'ascii',overwrite=True)

