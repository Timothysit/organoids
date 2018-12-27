import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 
import pandas as pd
from matplotlib import rcParams


ratePreTTX = np.array([0.02, 0.02, 0.02, 0.02, 0.32, 0.02, 0.02, 0.02])
ratePostTTX = np.array([0, 0, 0, 0, 0, 0, 0, 0])
spikeRate = np.concatenate((ratePreTTX, ratePostTTX))
condition = np.array(['pre-TTX', 'post-TTX'])
condition = np.repeat(condition, [8, 8], axis = 0)

# ratePreTTX_sem = np.std(ratePreTTX) / np.sqrt(len(ratePreTTX))
# ratePostTTX_sem = np.std(ratePostTTX) / np.sqrt(len(ratePostTTX))
# sem = [ratePreTTX_sem, ratePostTTX_sem]
## Supp figure 2l 

df = {'Condition': condition, 
'Spike Frequency': spikeRate}

df = pd.DataFrame(data = df)

# plt.figure(figsize=(4,10))
# ax=plt.subplots(111)
# figure size in inches
sns.set(rc={'figure.figsize':(6, 6 * 1.6)})
sns.set_style("white")
sns.set_style("ticks")
sns.barplot(x='Condition', y='Spike Frequency', ci=68, data=df)
# ci=68 uses the SEM instead of STD
# sns.boxplot(x='Condition', y='Spike Frequency', data=df)
# sns.catplot(x='Condition', y='Spike Frequency', kind='bar', data=df)
# some reason, changing the figure size does not affect catplot
# sns.pointplot(x='Condition', y='Spike Frequency', data=df,
# 	dodge=True, join=False)
sns.swarmplot(x='Condition', y='Spike Frequency', data=df, color='grey')


plt.ylabel('Spike Rate (Hz)')
plt.xlabel('')
sns.despine(top=True, right=True)
plt.savefig('supp_fig_2l_wTicks.eps', dpi=300)
plt.show()