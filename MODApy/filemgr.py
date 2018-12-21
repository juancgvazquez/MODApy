import pandas as pd
from os import path, remove
from MODApy.vcfmgr import ParsedVCF
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib_venn as venn

'''
Helper function to create Hyperlinks
'''


def make_hyperlink(value: str, urltype):
	if urltype == 'RSID':
		url = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs={}"
		if type(value) == str:
			return '=HYPERLINK("%s", "%s")' % (
				url.format(value.replace(';', ',').split(',')[0]), value.replace(';', ',').split(',')[0])
	if urltype == 'GENE':
		url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene={}"
		return '=HYPERLINK("%s", "%s")' % (url.format(value), value)
	else:
		return


'''
Helper function to get statisticts out of vcfs
'''


def getstats(self, type=0):
	stats = {}
	if (type == 0):
		if all(col in self.columns for col in ['ZIGOSITY', 'VARTYPE', 'PUTATIVE_IMPACT', 'EFFECT']):
			vcfstats = self.groupby([self.index.get_level_values(0), 'ZIGOSITY', 'VARTYPE', 'PUTATIVE_IMPACT',
									 'EFFECT']).size().to_frame(name='count')
			vcfstats.name = 'stats'
			stats['df'] = []
			stats['df'].append(vcfstats)
			plt.pie(vcfstats.groupby(level=0).count(), labels=vcfstats.groupby(level=0).count().index.values)
			my_circle = plt.Circle((0, 0), 0.7, color='white')
			chromVars = plt.gcf()
			chromVars.gca().add_artist(my_circle)
			stats['graphs'] = []
			stats['graphs'].append(chromVars)
	elif (type == 1):
		if 'TRIOS' in self.columns:
			trios = self.groupby('TRIOS', sort=False).size()
			A, B, C = self.name.split(':')
			triosgraph = plt.figure()
			venn.venn3(trios, set_labels=[A, B, C], set_colors=['b', 'r', 'g'])
			triosgraph.savefig('./triostmp.png', dpi=triosgraph.dpi)
			triosgraph.clf()
		elif 'DUOS' in self.columns:
			duos = self.groupby('DUOS', sort=False).size()
			A, B = self.name.split(':')
			duosgraph = plt.figure()
			venn.venn2(duos, set_labels=[A, B], set_colors=['b', 'r'])
			duosgraph.savefig('./duostmp.png', dpi=duosgraph.dpi)
			duosgraph.clf()
	return stats


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


def df_to_excel(df1: ParsedVCF, outpath):
	if (len(df1.index) > 65300):
		output = pd.ExcelWriter(outpath, engine='xlsxwriter', options={'strings_to_urls': False})
	else:
		output = pd.ExcelWriter(outpath)
	print('changing ID to url')

	try:
		df1['ID'] = df1['ID'].apply(lambda x: make_hyperlink(x, 'RSID'))
		df1['GENE_ID'] = df1['GENE_ID'].apply(lambda x: make_hyperlink(x, 'GENE'))
	except:
		print('Cant parse ID Field')

	# temp column drop until applied in config
	df1.drop(columns=['DISTANCE', 'GENE_NAME', 'ERRORS / WARNINGS / INFO'], inplace=True)
	# reordering columns so ID and GENE ID are first
	cols_selected = ['GENE_ID', 'ID', 'HGVS.P', 'HGVS.C', 'EFFECT', 'PUTATIVE_IMPACT', 'VARTYPE', '1000GP3_AF',
					 '1000GP3_AFR_AF', '1000GP3_AMR_AF', '1000GP3_EAS_AF', '1000GP3_EUR_AF', '1000GP3_SAS_AF',
					 'ESP6500_MAF_EA', 'ESP6500_MAF_AA', 'ESP6500_MAF_ALL', 'CLINVAR_CLNSIG', 'CLINVAR_CLNDSDB',
					 'CLINVAR_CLNDSDBID', 'CLINVAR_CLNDBN', 'CLINVAR_CLNREVSTAT', 'CLINVAR_CLNACC', 'ZIGOSITY',
					 'PolyPhen_Pred', 'PolyPhen_Score']
	# collist = [x for x in df1.columns if x not in cols_selected]
	df1.sort_index(inplace=True)
	df1 = df1[cols_selected]
	# just for fleni, temp until i do it form cfg

	# removing qual column if duos or trios
	# singlecols = ['QUAL']
	# if df1.columns[-1] == 'TRIOS' or df1.columns[-1] == 'DUOS':
	#	df1.drop(columns=singlecols, inplace=True)
	workbook = output.book
	datasheet = workbook.add_worksheet('DATA')
	statsheet = workbook.add_worksheet('STATISTICS')
	output.sheets['DATA'] = datasheet
	format1 = workbook.add_format({'num_format': '###,###,###'})
	hyper = workbook.add_format({'hyperlink': True})
	datasheet.set_column('B:B', 18, format1)
	datasheet.set_column('E:F', 18, hyper)
	df1.to_excel(output, sheet_name='DATA', merge_cells=False)
	stats = getstats(df1)
	output.sheets['STATISTICS'] = statsheet
	for i in range(len(stats['df'])):
		stats['df'][i].to_excel(output, sheet_name='STATISTICS')
	for i in range(len(stats['graphs'])):
		stats['graphs'][i].savefig('./tempgraph.png')
		statsheet.insert_image('H2', './tempgraph.png')
	if path.isfile('./duostmp.png'):
		statsheet.insert_image('H25', './duostmp.png')
	elif path.isfile('./triostmp.png'):
		statsheet.insert_image('H25', './triostmp.png')
	output.save()
	try:
		remove('./tempgraph.png')
	except:
		print('couldnt remove tempgraph.png')
	try:
		remove('./triostmp.png')
	except:
		print('couldnt remove triostmp.png')

	try:
		remove('./duostmp.png')
	except:
		print('couldnt remove duostmp.png')


'''
Saving new VCF (WIP)
'''


def df_to_vcf(df1: pd.DataFrame, outpath: str):
	header = """##fileformat=VCFv4.1
##source=MODApy
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
