import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
data = [ [0]*3 for i in range(3)]
with open('schedulerdata.txt') as f:
    for line in f: # read rest of lines
        data.append([int(x) for x in line.split()])
x = np.array([row[0] for row in data])
y1 = np.array([row[1] for row in data])
y2 = np.array([row[2] for row in data])
plt.rcParams["font.family"] = "sans"
fig = plt.figure()
ax = plt.gca()
ax.margins(x=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.scatter(x,y1,color='xkcd:sky blue',alpha=0.3,edgecolors='none', label='ConvScheduler', zorder=1)
ax.scatter(x,y2,color='xkcd:light red',alpha=0.3,edgecolors='none', label='SolvScheduler', zorder=1)
ax.set_yscale('log')
ax.set_ylim(top=1e6)
plt.xticks(np.arange(0, 1001, step=100))
def onepoint66(x, a):
    return a * np.power(x, 5/3)
def onepoint4(x, a):
    return a * np.power(x, 7/5)
newX = np.linspace(0, 1000, 1001)
popt, pcov = curve_fit(onepoint66, x, y1)
plt.plot(newX, onepoint66(newX, popt[0]), color='xkcd:darkish blue', alpha=0.7, lw=3, zorder=2, label="{0:.3f}x^(5/3), R^2 = {1:.3f}".format(popt[0], r2_score(y1, onepoint66(x, popt[0]))))
popt, pcov = curve_fit(onepoint4, x, y2)
plt.plot(newX, onepoint4(newX, popt[0]), color='xkcd:darkish red', alpha=0.7, lw=3, zorder=2, label="{0:.3f}x^(7/5), R^2 = {1:.3f}".format(popt[0], r2_score(y2, onepoint4(x, popt[0]))))
leg = plt.legend(loc='lower right')
for lh in leg.legend_handles: 
    lh.set_alpha(0.8)
plt.xlabel("P (sum of processing times of jobs)", labelpad=15)
plt.ylabel("Average Runtime (μs)", labelpad=15)
plt.title("Runtime of Scheduling Algorithm Implementations", fontweight="bold", pad=20)
plt.savefig('figure.png', dpi=300, bbox_inches="tight")
plt.show()