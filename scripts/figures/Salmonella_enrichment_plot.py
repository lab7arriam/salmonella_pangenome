import pandas as pd
import string
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab


sns.set(style="ticks")
sns.set_style("darkgrid")

params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}

pylab.rcParams.update(params)

# read in the table as a pandas dataframe
df_go = pd.read_table("./GO_topGO_enrichments_main.csv", sep=',', header=0)

df_go['Significant'] = df_go['Significant'].astype(int)
df_go['Annotated'] = df_go['Annotated'].astype(int)

df_go['division'] = df_go['Significant'] / df_go['Annotated']


adj_pic_vec = ['Bird_adj_pos', 'Cattle_adj_pos', 'Pig_adj_pos', 'Human_adj_pos', 'Virulence']
raw_pic_vec = ['Bird_raw_neg', 'Cattle_raw_neg', 'Pig_raw_neg', 'Human_raw_neg']
raw_pic_vec_pos = ['Bird_raw_pos', 'Cattle_raw_pos', 'Pig_raw_pos', 'Human_raw_pos', 'Virulence']

# Filter the data by the vector of group names
GO_results_adj = df_go[df_go['group'].isin(adj_pic_vec)]

GO_results_raw = df_go[df_go['group'].isin(raw_pic_vec)]

GO_results_raw_pos = df_go[df_go['group'].isin(raw_pic_vec_pos)]

# Select the top 7 enriched terms per group and ontology
Top_enrich_sig_adj = (GO_results_adj.sort_values(by=['Significant', 'classicFisher'], ascending=False)
                      .groupby(['group', 'ontology'])
                      .head(7))

Top_enrich_sig_adj['group'] = Top_enrich_sig_adj.group.apply(lambda x: ''.join(x.split('_')[0]))


# subset the dataframe by ontology
bp = Top_enrich_sig_adj[Top_enrich_sig_adj['ontology'] == 'GO:biological_process']
mf = Top_enrich_sig_adj[Top_enrich_sig_adj['ontology'] == 'GO:molecular_function']
cc = Top_enrich_sig_adj[Top_enrich_sig_adj['ontology'] == 'GO:cellular_component']


# create the figure and axes
fig, axs = plt.subplots(ncols=3, figsize=(25,18))
# fig, axs = plt.subplots(figsize=(4,13))


# plot each ontology as a scatterplot
sns.scatterplot(data=bp, x='group', y='Term', size='division', legend=False, hue = 'classicFisher', sizes=(100,250), ax=axs[0], palette='Blues_r').set(xlabel='', ylabel='')
sns.scatterplot(data=mf, x='group', y='Term', size='division', legend=False, hue = 'classicFisher', sizes=(100,250), ax=axs[1], palette='Blues_r').set(xlabel='', ylabel='');
sns.scatterplot(data=cc, x='group', y='Term', size='division', legend=True, hue = 'classicFisher', ax=axs[2], sizes=(100,250), palette='Blues_r').set(xlabel='', ylabel='');


handles, labels = axs[2].get_legend_handles_labels()

labels[0] = 'p-value'
labels[6] = 'Enrichment ratio'
fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1))


for n, ax in enumerate(axs):
    ax.text(0, 1.005, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')
    ax.tick_params('x',labelrotation=50)
    
plt.legend([],[], frameon=False)


fig.text(0.5, 0, 'Group', ha='center', fontsize=30)

# set common y-axis label
fig.text(0, 0.5, 'Term', va='center', rotation='vertical', fontsize=30)

plt.savefig('picture.png', bbox_inches='tight')

fig.tight_layout()
