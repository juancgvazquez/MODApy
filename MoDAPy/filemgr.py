import pandas as pd
from os import path

'''
Helper function to create Hyperlinks
'''
def make_hyperlink(value,urltype):
	if urltype == 'RSID':
		url = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs={}"
		return '=HYPERLINK("%s", "%s")' % (url.format(value), value)
	elif urltype == 'GENE':
		url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene={}"
		return '=HYPERLINK("%s", "%s")' % (url.format(value), value)
	else:
		return



'''
VCF DataFrame to Xlsx
'''

def checkFile(filePath, extension):
	if path.isfile(filePath):
		fileName, fileExtension = path.splitext(filePath)
		if extension == fileExtension:
			return True

	print(filePath, "couldn't be found. Please check if file exists and that it's extension is",
		  "'" + extension + "'")
	exit(1)


def df_to_excel(df1:pd.DataFrame, outpath):
	#output = pd.ExcelWriter(outpath)

	#	Aca convertiría el campo en enlace, pero todavía hay que evaluar que hacer con los múltiples rs
	if (len(df1.index) > 65300):
		output = pd.ExcelWriter(outpath, engine='xlsxwriter', options={'strings_to_urls': False})
	else:
		output = pd.ExcelWriter(outpath)
	print('changing ID to url')

	try:
		df1['ID'] = df1['ID'].apply(lambda x: make_hyperlink(x,'RSID'))
		df1['Gene_Name'] = df1['Gene_Name'].apply(lambda x: make_hyperlink(x,'GENE'))
	except:
		print('Cant parse ID Field')

	#temp column drop until applied in config
	df1.drop(columns=['Distance','Gene_ID', 'ERRORS / WARNINGS / INFO'])
	df1.sort_index(inplace=True)
	df1.to_excel(output, sheet_name='Result')
	workbook = output.book
	worksheet = output.sheets['Result']
	format1 = workbook.add_format({'num_format': '###,###,###'})
	worksheet.set_column('B:B', 18,format1)
	output.save()


'''
Saving new VCF
'''


def df_to_vcf(df1: pd.DataFrame, outpath: str):
	header = """##fileformat=VCFv4.1
##source=MoDAPy
##reference=hg19
##INFO=<ID=ZIG,Number=.,Type=String,Description="Indicates Zigosity' ">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=MG_INDEL_Filter,Description="QD < 2.0  ||  FS > 200.0  ||  ReadPosRankSum < -20.0">
##FILTER=<ID=MG_SNP_Filter,Description="QD < 2.0  || FS > 60.0 || MQ < 40.0 ||  MQRankSum < -12.5  ||  ReadPosRankSum < -8.0">
#CHROM	POS	REF	ALT	ID	QUAL	FILTER	INFO
"""
	output_VCF = outpath
	with open(output_VCF, 'w') as cyvcf2:
		cyvcf2.write(header)

	info_list = df1.columns[4:20]
	df1['INFO'] = df1[info_list].apply(lambda x: 'ANN=' + '|'.join(x), axis=1)
	df1['INFO'] = 'ZIG=' + df1['ZIGOSITY'] + ';' + df1['INFO'] + ';TRI=' + df1['Trios']
	df1 = df1.drop(info_list, axis=1).drop(['ZIGOSITY', 'Trios'], axis=1)

	df1.reset_index().to_csv(output_VCF, sep='\t', mode='a', index=False)
