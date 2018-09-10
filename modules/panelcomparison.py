import pandas as pd
from sys import argv

def compare(file1, file2):
	print(file1)
	print(file2)
	df1 = pd.read_excel(file1)
	df2 = pd.read_excel(file2)

	set1 = df1['CHROM'].astype(str).str.cat(df1['POS'].values.astype(str), sep=':')
	set2 = df2['#CHROM'].astype(str).str.cat(df2['POS'].values.astype(str), sep=':')

	set1 = set(x.replace(',', '') for x in set1)
	set2 = set(x.replace(',', '') for x in set2)

	diferencia12 = [x for x in set1 if x not in set2]
	diferencia21 = [x for x in set2 if x not in set1]

	print("Datos del archivo 1 que no se encuentran en el 2:\n")
	print(diferencia12, '\n')
	print("Datos del archivo 2 que no se encuentran en el 1:\n")
	print(diferencia21)


if __name__ == '__main__':
	print("USAGE: First vcf with 'CHROM' second with '#CHROM'")
	compare(argv[1], argv[2])
