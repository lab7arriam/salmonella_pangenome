import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import string


# Set the color palette for all plots
colors = sns.color_palette("tab20", 9)
sns.set(style="ticks")
sns.set_style("darkgrid")

fig = plt.figure(figsize=(12, 10))

gs = fig.add_gridspec(2,2)

ax3 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax1 = fig.add_subplot(gs[1, :])


for n, ax in enumerate([ax3, ax2, ax1]):
    ax.text(-0.1, 1.1, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')


# Load data for the first plot
df_ratio = pd.read_csv('./host_ratio.csv')

# Create the violin plot using Seaborn for the first plot
sns.violinplot(x='host', y='ratio', data=df_ratio, palette=[colors[1], colors[0], colors[3], colors[2]], ax=ax1)

# Add titles and labels for the first plot
ax1.set_xlabel('Host', fontsize=12)
ax1.set_ylabel('Percentage of core genes', fontsize=12)

# Remove top and right spines for the first plot
sns.despine(ax=ax1)


# Load data for the second plot
roary = pd.read_table('./gene_presence_absence.csv', sep=',', low_memory=False)
roary.set_index('Gene', inplace=True)
roary.drop(list(roary.columns[:13]), axis=1, inplace=True)
roary.replace('.{2,100}', 1, regex=True, inplace=True)
roary.replace(np.nan, 0, regex=True, inplace=True)
idx = roary.sum(axis=1).sort_values(ascending=False).index
roary_sorted = roary.loc[idx]

# Create the histogram plot using Seaborn for the second plot
sns.histplot(data=roary.sum(axis=1), kde=False, bins=roary.shape[1],
             ax=ax2, color=colors[2], log_scale=(False,True))

# Set the logarithmic scale for y-axis and update the axis label for the second plot
ax2.set_yscale('log')
ax2.set_xlabel('Number of genomes', fontsize=12)
ax2.set_ylabel('Number of genes (log scale)', fontsize=12)

# Remove spines and add a grid for the second plot
sns.despine(ax=ax3, left=True, bottom=True)
ax3.grid(axis='y', alpha=0.3)


# Load data for the third plot
plotdf = pd.read_csv('./saturation_curve_data.csv')


# Create the line plot with fill between using Seaborn for the third plot
sns.lineplot(x='N', y='accessory size', data=plotdf, ax=ax3, color=colors[5])
ax3.fill_between(x=plotdf['N'], y1=plotdf['accessory size'] - plotdf['std'],
                 y2=plotdf['accessory size'] + plotdf['std'], color=colors[5], alpha=0.5)

# Add titles and labels for the third plot
ax3.set_xlabel('Number of genomes', fontsize=12)
ax3.set_ylabel('Number of genes', fontsize=12)

# Remove top and right spines for the third plot
sns.despine(ax=ax3)


# Adjust the spacing between the subplots
plt.subplots_adjust(hspace=0.3, wspace=0.3)


plt.savefig('pangenome_properties.png', dpi=512)


# Show the plot
plt.show()
