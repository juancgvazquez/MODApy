from collections import OrderedDict
import cyvcf2
import pandas as pd

'''
Simple parser, Uses cyvcf2 Returns a Reader Object
'''


def simple_parser(vcf):
	try:
		pvcf = cyvcf2.Reader(vcf)
	except:
		print('error loading vcf')

	return pvcf


'''
Detailed parser, Uses cyvcf2 and returns a DataFrame containing all info.
'''


def vcf_to_df(vcf):
	try:
		pvcf = cyvcf2.Reader(vcf)
	except:
		print('error loading vcf')


	info_dict = OrderedDict()
	headers = {'CHROM': [], 'POS': [], 'REF': [], 'ALT': [], 'ID': [], 'QUAL': [], 'FILTER': [], 'ZIGOSITY': []}
	vcf_dict = OrderedDict(headers.copy())
	infokeylist = ['AC', 'SAMPLES_AF', 'AN', 'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'QD', 'SOR', 'dbSNPBuildID', 'ANN']
	milgp3keylist = ['1000Gp3_AF','1000Gp3_AFR_AF','1000Gp3_AMR_AF','1000Gp3_EAS_AF','1000Gp3_EUR_AF','1000Gp3_SAS_AF']
	esp6500keylist = ['ESP6500_MAF_EA','ESP6500_MAF_AA','ESP6500_MAF_ALL','ESP6500_PH']
	clinvarkeylist = ['CLINVAR_CLNSIG','CLINVAR_CLNDSDB','CLINVAR_CLNDSDBID','CLINVAR_CLNDBN','CLINVAR_CLNREVSTAT','CLINVAR_CLNACC']


	for variant in pvcf:
		for (key, value) in variant.INFO:
			if key in infokeylist:
				if key in info_dict:
					info_dict[key].append(value)
				else:
					info_dict[key] = [value]
			elif key in milgp3keylist:
				if key in info_dict:
					info_dict[key].append(value)
				else:
					info_dict[key] = [value]
			elif key in esp6500keylist:
				if key in info_dict:
					info_dict[key].append(value)
				else:
					info_dict[key] = [value]
			elif key in clinvarkeylist:
				if key in info_dict:
					info_dict[key].append(value)
				else:
					info_dict[key] = [value]

		zigosity = ''
		if variant.gt_types == 0:
			zigosity = 'HOM_REF'
		elif variant.gt_types == 1:
			zigosity = 'HET'
		elif variant.gt_types == 2:
			zigosity = 'UNKNOWN'
		elif variant.gt_types == 3:
			zigosity = 'HOM_ALT'
		for header in vcf_dict:
			if header != 'ZIGOSITY':
				vcf_dict[header].append(getattr(variant, header))
			else:
				vcf_dict['ZIGOSITY'].append(zigosity)

	maindf = pd.DataFrame.from_dict(data=vcf_dict, orient='index').transpose()
	try:
		maindf['ALT'] = maindf['ALT'].str.join(',')
	except:
		pass
	maindf.index.set_names('Variant', inplace=True)
	if info_dict['ANN']:
		annlist = info_dict.pop('ANN')
		# Preparo lista de cabeceras
		annhead = pvcf.get_header_type('ANN')['Description'].strip('"Functional annotations:"')
		annheaderlist = [x.strip() for x in annhead.split('|')]
		annlsplit = [x.split(',') for x in annlist]
		anndf = pd.DataFrame(annlsplit).stack()
		anndf = anndf.str.split('|', expand=True)
		anndf.columns = annheaderlist
		anndf.index.set_names(['Variant', 'Ann_N°'], inplace=True)
	else:
		print('no annotations')
		anndf = pd.DataFrame({'Variant': []})

	infodf = pd.DataFrame.from_dict(data=info_dict, orient='index').transpose()
	infodf.index.set_names('Variant', inplace=True)
	fulldf = maindf.join(infodf, how='inner').join(anndf, how='inner')
	fulldf.reset_index(inplace=True)
	fulldf = fulldf.drop(['Variant', 'Ann_N°'], axis=1)
	fulldf.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True)
	fulldf.columns = fulldf.columns.str.strip("'")
	try:
		fulldf.name = pvcf.samples[0]
	except:
		fulldf.name = ''
	return fulldf


