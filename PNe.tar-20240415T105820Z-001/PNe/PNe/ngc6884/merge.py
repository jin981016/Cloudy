import numpy as np
import pandas as pd

d1 = pd.read_csv('cngc6884_3600s.0054.txt', sep='  ',names = ['wavelength','flux'])
d2 = pd.read_csv('cngc6884_3600s.0055.txt', sep='  ',names = ['wavelength','flux'])
dd1 = pd.DataFrame(d1,columns = ['wavelength','flux'])
dd2 = pd.DataFrame(d2,columns = ['wavelength','flux'])

df = pd.concat([dd1,dd2],axis=0)
df = df.sort_values(by=['wavelength'])
df.to_csv('alpha.txt', header = None)