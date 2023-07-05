import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('./assemblies_clusters_300.csv')

df['GC'] = df['GC'].round(1)

df = df[df['Hypothetical Proteins'] < 4000]


sns.set(style="ticks")
sns.set_style("darkgrid")

plt.figure(figsize=(12, 10))


ax = sns.scatterplot(data=df, x='Genome_length', y="Hypothetical Proteins", hue='GC', palette="deep", s=100, legend="brief")

sns.regplot(data=df, x='Genome_length', y="Hypothetical Proteins", scatter=False, ax=ax, ci=95)

ax.set_xlabel('Genome length', fontsize=20)
ax.set_ylabel('Hypothetical Proteins', fontsize=20)

ax.legend(title="GC", fontsize=15, bbox_to_anchor=[1.15, 1.014])

plt.savefig('length_proteins_gc.png', dpi=512)
plt.show()
