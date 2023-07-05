import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

colors = sns.color_palette("tab20", 9)
sns.set(style="ticks")
sns.set_style("darkgrid")

# read in the table as a pandas dataframe
df_go = pd.read_table("./GO_topGO_enrichments_main.csv", header=0)


df_go['Significant'] = df_go['Significant'].astype(int)
df_go['Annotated'] = df_go['Annotated'].astype(int)

df_go['division'] = df_go['Significant'] / df_go['Annotated']


adj_pic_vec = ['Bird_adj_pos', 'Cattle_adj_pos', 'Pig_adj_pos', 'Human_adj_pos'] #, 'Virulence']
raw_pic_vec = ['Bird_raw_neg', 'Cattle_raw_neg', 'Pig_raw_neg', 'Human_raw_neg']
raw_pic_vec_pos = ['Bird_raw_pos', 'Cattle_raw_pos', 'Pig_raw_pos', 'Human_raw_pos', 'Virulence']

# Filter the data by the vector of group names
GO_results_adj = df_go[df_go['group'].isin(adj_pic_vec)]

GO_results_raw = GO_results_all[GO_results_all['group'].isin(raw_pic_vec)]
GO_results_raw_pos = GO_results_all[GO_results_all['group'].isin(raw_pic_vec_pos)]

# Select the top 7 enriched terms per group and ontology
Top_enrich_sig_adj = (GO_results_adj.sort_values(by=['Significant', 'classicFisher'], ascending=False)
                      .groupby(['group', 'ontology'])
                      .head(7))

#Top_enrich_sig_adj['group'] = Top_enrich_sig_adj.group.apply(lambda x: ''.join(x.split('_')[0]))

# subset the dataframe by ontology
bp = Top_enrich_sig_adj[Top_enrich_sig_adj['ontology'] == 'GO:biological_process']
mf = Top_enrich_sig_adj[Top_enrich_sig_adj['ontology'] == 'GO:molecular_function']
cc = Top_enrich_sig_adj[Top_enrich_sig_adj['ontology'] == 'GO:cellular_component']

# create the figure and axes
fig, axs = plt.subplots(ncols=3, figsize=(25,18))

for n, ax in enumerate(axs):
    ax.text(0, 1.005, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)

# plot each ontology as a scatterplot
sns.scatterplot(data=bp, x='group', y='Term', size='division', sizes=(50,500), ax=axs[0], legend=False).set(xlabel='', ylabel='')
sns.scatterplot(data=mf, x='group', y='Term', size='division', sizes=(50,500), ax=axs[1], legend=False, color='#33a02c').set(xlabel='', ylabel='');
sns.scatterplot(data=cc, x='group', y='Term', size='division', sizes=(50,500), ax=axs[2], legend=False, color='#ff7f0e').set(xlabel='', ylabel='');

# set the titles for each axis
# axs[0].set_title('Biological Process')
# axs[1].set_title('Molecular Function')
# axs[2].set_title('Cellular Component')

fig.text(0.5, 0.08, 'Group', ha='center', fontsize='large')

# set common y-axis label
fig.text(-0.04, 0.5, 'Term', va='center', rotation='vertical', fontsize='large')


# # adjust the spacing between subplots
plt.subplots_adjust(wspace=0.9)


# display the figure
plt.show()
