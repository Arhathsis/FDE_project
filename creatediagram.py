import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

table = pd.read_csv('file1.txt', header=None, sep=' ')
table2 = np.array(2 - table)
#med = np.median(table2[np.where(~np.isnan(table2))])
#table3 = np.copy(table2)
#table3[np.where(np.isnan(table3))] = med


x = np.array([2000, 4472, 10000, 22361, 50000])
y = np.array([5, 9, 16, 28, 50])

xx, yy = np.meshgrid(x, y)

fig, ax = plt.subplots()
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel(r'$N$')
ax.set_ylabel(r'$l_s$')
color_map = plt.cm.get_cmap('Blues')
reversed_color_map = color_map.reversed()
cmesh1 = ax.pcolormesh(np.log10(xx), np.log10(yy), table2, shading='gouraud', cmap=reversed_color_map)
cbar = fig.colorbar(cmesh1, ax=ax)
cbar.set_label(r'$\Delta D$', rotation=270)
ax.set_xticks(ticks=np.log10(x))
ax.set_xticklabels(x)
ax.set_yticks(ticks=np.log10(y))
ax.set_yticklabels(y)
ax.set_title('FDE 2d dispersion distribution')
fig.tight_layout()
