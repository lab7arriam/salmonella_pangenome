import numpy as np
import pandas as pd
from collections import Counter

from sklearn.preprocessing import MultiLabelBinarizer

groups = {
    'Bird': ['Anas platyrhynchos', 'Avian', 'Chroicocephalus novaehollandiae',
             'Cormorant (unspecified)', 'Finch', 'Gallus gallus', 'Goose',
             'Gull (unspecified)', 'Hawk', 'Meleagris gallopavo', 'Parrot',
             'Passer domesticus', 'Passer hispaniolensis', 'Zonotrichia leucophrys'],
    'Equine': ['Equus africanus asinus', 'Equus ferus caballus'], #deleted
    'Sheep': ['Ovis aries'], #deleted
    'Cattle': ['Bos grunniens', 'Bos taurus'],
    'Pig': ['Sus scrofa domesticus'],
    'Goat': ['Capra aegagrus'], #deleted
    'Canidae': ['Canis latrans', 'Canis lupus familiaris'], #deleted
    'Reptiles': ['Chamaeleonidae', 'Crocodylus', 'Iguana',
                 'Lizard (Cordylus niger)', 'Lacertilia', 'Pogona vitticeps',
                 'Reptile (unspecified)'], #deleted 
    'Environment': ['Lab', 'environmental', 'Fish', 'Crustacean'], #deleted
    'Cat': ['Felis catus domesticus'], #deleted
    'Human': ['Homo sapiens'],
    'Rodent': ['Mus musculus'], #deleted
    'Marsupialia': ['Opossum'] #deleted
}

df_ncbi = pd.read_csv('df_biosample.csv')
serovars = pd.read_csv('ncbi_serovar.csv')
df_host = df_ncbi[['host', 'NCBI', 'serovar', 'isolation_source']]

order = serovars.serovars.to_list()

former_dict = {}

for group in list_of_assembly_dereplicated_groups:
    group_hosts = []

    for item in group:
        item = '.'.join((item.split('.'))[0:2])
        host = df_host.loc[df_host['NCBI'] == item, 'host'].iloc[0]
        if pd.isna(host):
            host = 'missing'
        group_hosts.append(host)
    
    former_dict.update({'.'.join((group[0].split('.'))[0:2]) : list(set(group_hosts))})
    
    
new_dict = {i: former_dict.get(i) for i in order}
all_val = [(k, v) for k in former_dict for v in former_dict[k]]

        
df = pd.DataFrame(all_val, columns=['key','val']).set_index('key')
df_count = df.pivot_table(index='key', columns='val', aggfunc=len)

mlb = MultiLabelBinarizer()
df = pd.DataFrame(mlb.fit_transform(former_dict.values()),columns=mlb.classes_, index=former_dict.keys())
