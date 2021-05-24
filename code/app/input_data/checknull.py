import pandas as pd

df = pd.read_csv('clean_yeast_data.csv')

print(df.isnull().sum())


