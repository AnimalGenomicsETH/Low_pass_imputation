## FIGURE 1A

# Load packages
from upsetplot import UpSet, from_memberships
import matplotlib as mpl

# Load data
indicators = [['DV','GATK'],['DV'],['GATK']],[15730925,1709313,1470358]
upset = from_memberships(*indicators)

# Divide by 1 million so sizes because million of SNPs
upset/=1e6

# pick scale for colorbar and colormap name
norm = mpl.colors.Normalize(vmin=1,vmax=3)
cmap = 'cividis'

# make the figure, add the upsetplot
fig = plt.figure(figsize=(9, 6))
g = UpSet(a,sort_by='cardinality',show_counts='%.2f',element_size=None)

# colour the specific bars by the *HARDCODED* Ti:Tv ratio
g.style_subsets(present="DV", absent="GATK", facecolor=mpl.cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(2.07))
g.style_subsets(present=("DV",'GATK'), facecolor=mpl.cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(2.25))
g.style_subsets(present="GATK", absent="DV", facecolor=mpl.cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(1.44))

# neaten up plot and add colorbar
axes = g.plot(fig=fig)
axes['totals'].axis('off')
axes['intersections'].set_ylabel('biSNPs (millions)')
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes['intersections'],orientation='horizontal',location='top',pad=0.1,label='Ti:Tv')

fig.savefig('figure_1a.svg')
