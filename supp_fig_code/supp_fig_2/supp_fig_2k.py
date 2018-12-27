import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 
import pandas as pd
from matplotlib import rcParams


condition = np.array(['pre-TTX', 'post-TTX'])
activeElectrodes = np.array([8, 0])


## Supp figure 2l 

df = {'Condition': condition, 
'electrodeCount': activeElectrodes}

df = pd.DataFrame(data = df)

# plt.figure(figsize=(4,10))
# ax=plt.subplots(111)
# figure size in inches
sns.set(rc={'figure.figsize':(6, 8)})
sns.set_style("white")
sns.set_style("ticks")
sns.barplot(x='Condition', y='electrodeCount', data=df)
plt.ylabel('Number of active electrodes')
plt.xlabel('')
sns.despine(top=True, right=True)
plt.savefig('supp_fig_2k_wTicks.eps', dpi=300)
plt.show()