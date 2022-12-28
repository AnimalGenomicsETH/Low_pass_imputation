

# Load data
    
indicators = [['Chip','DV','GATK'],['DV','Chip'],['GATK','Chip'],['Chip']],[492229,959,79,5309]
upset = from_memberships(*indicators)

fig = mpl.pyplot.figure(figsize=(9, 6))
axes = UpSet(upset,sort_by='cardinality',facecolor='gray',show_counts='%d',element_size=None).plot(fig=fig)
axes['intersections'].set_ylabel('biSNPs')
axes['intersections'].set_yscale('log')

HIGH = '#F10E3C'
MODERATE = '#EC8013'
LOW = '#D2CD2D'

axes['intersections'].text(0,1e5,'4461',zorder=10,ha='left',c=HIGH)
axes['intersections'].text(0,3e4,'1855',zorder=10,ha='left',c=MODERATE)
axes['intersections'].text(0,9e3,'104',zorder=10,ha='left',c=LOW)

axes['intersections'].text(2,1e3,'5',zorder=10,ha='left',c=HIGH)
axes['intersections'].text(2,3e2,'2',zorder=10,ha='left',c=MODERATE)
axes['intersections'].text(2,9e1,'4',zorder=10,ha='left',c=LOW)

axes['intersections'].text(3,1e3,'0',zorder=10,ha='left',c=HIGH)
axes['intersections'].text(3,3e2,'1',zorder=10,ha='left',c=MODERATE)
axes['intersections'].text(3,9e1,'1',zorder=10,ha='left',c=LOW)

axes['intersections'].text(1,1e4,'HIGH',zorder=10,ha='left',c=HIGH)
axes['intersections'].text(1,3e3,'MODERATE',zorder=10,ha='left',c=MODERATE)
axes['intersections'].text(1,9e2,'LOW',zorder=10,ha='left',c=LOW)

axes['intersections'].text(1,9e4,'/',zorder=10,ha='left',c='k')

axes['totals'].axis('off')

fig.savefig('figure_1b.svg')
