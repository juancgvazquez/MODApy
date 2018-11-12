from collections import OrderedDict
import cyvcf2
import pandas as pd


class ParsedVCF(pd.DataFrame):
	def getstats(self):
		vcfdfstats = self.groupby(['ZIGOSITY', 'VARTYPE', 'Annotation_Impact', 'Annotation']).size().to_frame(
			name='count')

	@classmethod
	def from_vcf(cls, vcf):
		try:
			pVCF = cyvcf2.Reader(vcf)
		except:
			print('error loading vcf')
			exit(1)

		variantsDict = OrderedDict()
		for variant in pVCF:
			variantsDict[variant.CHROM + '_' + str(variant.POS) + '_' + variant.REF + '_' + ','.join(variant.ALT)] = {
				'ID': variant.ID, 'QUAL': variant.QUAL, 'FILTER': variant.FILTER}
			variantsDict[
				variant.CHROM + '_' + str(variant.POS) + '_' + variant.REF + '_' + ','.join(variant.ALT)].update(
				{k: v for (k, v) in variant.INFO})

		df1 = pd.DataFrame.from_dict(variantsDict, orient='index')
		if 'ANN' in df1.columns:
			anndf = df1['ANN']
			annhead = pVCF.get_header_type('ANN')['Description'].strip('"Functional annotations: \'"')
			annheaderlist = [x.strip() for x in annhead.split('|')]
			anndf = anndf.str.split(',', expand=True).stack()
			anndf = anndf.str.split('|', expand=True)
			anndf.columns = annheaderlist
			df1.drop(columns=['ANN'], inplace=True)
			anndf.index = anndf.index.droplevel(1)
			vcfdf = df1.join(anndf, how='inner')
		else:
			vcfdf = df1

		vcfdf.index = vcfdf.index.str.split('_', expand=True)
		vcfdf.index.names = ['CHROM', 'POS', 'REF', 'ALT']
		if 'HOM' in vcfdf.columns:
			vcfdf['HOM'] = vcfdf['HOM'].map({True: 'HOM'})
			vcfdf.HOM.fillna('HET', inplace=True)
			vcfdf.rename(columns={'HOM': 'ZIGOSITY'}, inplace=True)
		try:
			vcfdf.name = pVCF.samples[0]
		except:
			vcfdf.name = vcf.split('/')[-1]

		return vcfdf
