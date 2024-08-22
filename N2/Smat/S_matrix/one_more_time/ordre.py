import pandas as pd
from functools import reduce

my_folder_names = ['01','02','03','04','05']
geom = []

for folder in my_folder_names:
	my_columns = ['1', '2', '3', '4', '5', '6', '7', '8', 'val', '10', '11', 'state no.', '-et', 'symetrie']
	df = pd.read_csv('./'+folder+'/results/target.prop', names = my_columns, sep='\s{1,}', engine='python')
	int_df = df.drop(columns=['1', '2', '3', '4', '5', '6', '7', '8', '10', '11', 'state no.', '-et'])
	int_2_df = int_df[int_df.symetrie.notnull()]
	int_3_df = int_2_df.sort_values(by=['symetrie', 'val'])
	rslt_df = int_3_df.drop_duplicates('symetrie', keep='first', ignore_index=True)
	pd.set_option('display.max_rows', None)

	rslt_df.to_csv(folder+'_target_NRJ.csv', index=False)

	with open('./'+folder+'/Q.txt') as f:
		geom.append(f.readlines())
		geom[-1] = geom[-1][0][:-2]

df1 = pd.read_csv('01_target_NRJ.csv', engine='python')
df2 = pd.read_csv('02_target_NRJ.csv', engine='python')
df3 = pd.read_csv('03_target_NRJ.csv', engine='python')
df4 = pd.read_csv('04_target_NRJ.csv', engine='python')
df5 = pd.read_csv('05_target_NRJ.csv', engine='python')
dfs = [df1,df2,df3,df4,df5]

result = pd.merge(left=df1,right=df2,on=['symetrie'], how='outer', suffixes=('1','2'))
result = pd.merge(left=result,right=df3,on=['symetrie'], how='outer', suffixes=(None,'3'))
result = pd.merge(left=result,right=df4,on=['symetrie'], how='outer', suffixes=(None,'4'))
result = pd.merge(left=result,right=df5,on=['symetrie'], how='outer', suffixes=(None,'5'))

result = result.drop(columns=['symetrie'])

min = min(result.min(axis='columns'))

result = result - min

result = result.multiply(27.2114)

result = result + 15.43158 #minimum of N2+_pot_001.dat

result = pd.DataFrame(result.values)

result.columns = geom

result = result.T

result.to_csv('target_NRJ.csv', header=False, sep=' ')
