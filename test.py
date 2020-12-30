import pandas as pd
import numpy as np


d = {'col1': [1, 4, 3, 4, 5], 'col2': [4, 5, 6, 7, 8], 'col3': [7, 8, 9, 0, 1]}
df = pd.DataFrame(data=d)
print("Original DataFrame")
print(df)
print('Data from new_file.csv file:')
df.to_csv('new_file.csv', sep='\t', index=False)
new_df = pd.read_csv('new_file.csv')
print(new_df)