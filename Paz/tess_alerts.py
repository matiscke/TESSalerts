from astropy.table import Table
import numpy as np
from astropy.table import Column
from astropy import units as u   
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
import pandas as pd
import astroplan
import datetime as dt
from astroplan import Observer 
from astroplan import FixedTarget
import requests
from lxml import html

print("insert username:")
user = input('')
print("insert password:")
password = input('')

def download_tess_alert_table(user,password,table_name = 'alerts_tess_table.dat'):
    LOGIN_URL = 'https://tev.mit.edu/user/login/?next=http%3A//127.0.0.1%3A3032/toi/alerts/all/all/csv/'
    URL = 'https://tev.mit.edu/toi/alerts/all/all/csv/'
    session_requests = requests.session()
    result = session_requests.get(LOGIN_URL)
    tree = html.fromstring(result.text)
    authenticity_token = list(set(tree.xpath("//input[@name='csrfmiddlewaretoken']/@value")))[0]
    payload = {
        "username": user, 
        "password": password, 
        "csrfmiddlewaretoken": authenticity_token
              }
    result = session_requests.post(LOGIN_URL, data = payload, headers = dict(referer = LOGIN_URL))
    
    result = session_requests.get(URL, headers = dict(referer = URL))
    f = open(table_name,'w')
    f.write(str(result.content).replace(',','\t').replace(' ','_').replace('\\r\\n','\n').replace('\t\t','\t-\t').replace('\t\t','\t-\t').replace('\t\t','\t-\t').replace('\t\t','\t-\t').replace('\t\t','\t-\t')[2:-1]) 
    f.close()

download_tess_alert_table(user,password,table_name = 'alerts_tess_table.dat')

from astropy.utils.data import clear_download_cache 
clear_download_cache() 

#Read table with astropy
t_alerts =  Table.read('alerts_tess_table.dat',format = 'ascii')
alerts_ticid = t_alerts['tic_id'].data.data
alerts_toi_id = t_alerts['toi_id'].data.data
ra_alerts = t_alerts['RA']*u.deg
dec_alerts = t_alerts['Dec']*u.deg
coords_alerts = SkyCoord(ra_alerts,dec_alerts,unit = u.deg)
ra_alerts=coords_alerts.ra
dec_alerts=coords_alerts.dec
T_mag = t_alerts['Tmag'].data.data
Terr_mag = t_alerts['Tmag_Err'].data.data


#Check is the coordinates are visible from Calar Alto
lista_candidatas = []
indices_candidatas = []
for j,element in enumerate(coords_alerts):
    dec = element.dec
    if dec >-23*u.deg and dec <90*u.deg:
        lista_candidatas.append(element)
        indices_candidatas.append(j)



#Add visibility
months_list = []



#only candidates alerts (visible from Calar Alto)

T_mag_list=[]
Terr_mag_list=[]
RA_candidatas=[]
DEC_candidatas=[]
targets_candidatas=[]
toi_candidatas = []
p=len(lista_candidatas)
for i in range(p):
    print(str(i)+'/'+str(p))
    location = EarthLocation.from_geodetic(-2.546111*u.deg, 37.223611*u.deg, 2168*u.m) 
    obs = Observer(location=location, name="calar", timezone="CET")
    c=[astroplan.constraints.LocalTimeConstraint(min=dt.time(19,0), max = dt.time(5,0)),astroplan.constraints.AirmassConstraint(max=2) ] 
    indices= np.array(indices_candidatas) 
    targets_candidatas.append(alerts_ticid[indices[i]])
    toi_candidatas.append(alerts_toi_id[indices[i]])
    T_mag_list.append(T_mag[indices[i]])
    Terr_mag_list.append(Terr_mag[indices[i]])
    RA_candidatas.append(ra_alerts[indices[i]].value)
    DEC_candidatas.append(dec_alerts[indices[i]].value)
    coord =lista_candidatas[i]
    target = FixedTarget(coord=coord, name=alerts_ticid[indices[i]])
    months = astroplan.months_observable(constraints=c,observer=obs,targets=[target]) 
    months = str(months).replace(' ','').replace('[','').replace(']','') 
    months_list.append(months)

t_car_tess =  Table.read('colnames_carmenes_tess.dat',format = 'ascii')
name_t = t_car_tess['tess'].data
name_c = t_car_tess['Karmn']
c=np.array(targets_candidatas)
d=c.astype(str)
id_carmenes = []
for i,tic in enumerate(d):
    j = np.where(name_t==tic)[0]
    if len(j)>0:
        carmen_id = name_c[j[0]]
    else:
        carmen_id ='...'
    id_carmenes.append(carmen_id)

table=Table()
COL_targets=Column(data=targets_candidatas,name = 'tic_id',dtype = int)
table.add_column(COL_targets)
COL_toi=Column(data=toi_candidatas,name = 'toi_id',dtype = float)
table.add_column(COL_toi)
COl_ra=Column(data=RA_candidatas, name='RA', dtype= float)
table.add_column(COl_ra)
COl_dec=Column(data=DEC_candidatas, name='DEC', dtype= float)
table.add_column(COl_dec)
COl_Tmag=Column(data=T_mag_list, name='T_mag', dtype= float)
table.add_column(COl_Tmag)
COl_Terrmag=Column(data=Terr_mag_list, name='Tmag_Err', dtype= float)
table.add_column(COl_Terrmag)
COL_months = Column(data=months_list,name = 'Observable_Months',dtype = str)
table.add_column(COL_months)
COL_car = Column(data=id_carmenes,name = 'Karmn',dtype = str)
table.add_column(COL_car)
table.write('Alerts_TESS_Table.dat',format='ascii',overwrite=True)
