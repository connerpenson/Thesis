import pandas as pd
my_folder_names = ['01','02','03','04','05']
for folder in my_folder_names :
	my_columns = ['-et', 'symetrie', 'nrj']
	df = pd.read_csv('./'+folder+'/results/SelectedTargetStates.txt', header = None, names = my_columns, sep='\s{1,}', engine='python')
	pd.set_option('display.max_rows', None)
	rslt_df = df.sort_values(by=['symetrie', '-et', 'nrj'])
	print(rslt_df)

	rslt_df.to_csv(folder+'.csv', index=False)
