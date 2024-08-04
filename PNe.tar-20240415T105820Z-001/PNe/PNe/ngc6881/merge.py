import numpy as np
import pandas as pd

d1 = pd.read_csv('cngc6881_3300s.0026.txt', sep='  ',names = ['wavelength','flux'])
d2 = pd.read_csv('cngc6881_3300s.0029.txt', sep='  ',names = ['wavelength','flux'])
dd1 = pd.DataFrame(d1,columns = ['wavelength','flux'])
dd2 = pd.DataFrame(d2,columns = ['wavelength','flux'])

df = pd.concat([dd1,dd2],axis=0)
df.to_csv('beta.txt', header = None)