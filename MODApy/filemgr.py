import logging
from os import path, remove

import matplotlib
import pandas as pd

from MODApy.vcfmgr import ParsedVCF

matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib_venn as venn

logger = logging.getLogger(__name__)

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
        if all(col in self.columns for col in ['ZIGOSITY', 'VARTYPE', 'IMPACT', 'EFFECT']):
            vcfstats = self.groupby([self.index.get_level_values(0), 'ZIGOSITY', 'VARTYPE', 'IMPACT',
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

    logger.error(filePath, "couldn't be found. Please check if file exists and that it's extension is",
                 "'" + extension + "'")
    exit(1)


def df_to_excel(df1: ParsedVCF, outpath):
    if (len(df1.index) > 65300):
        output = pd.ExcelWriter(outpath, engine='xlsxwriter', options={'strings_to_urls': False})
    else:
        output = pd.ExcelWriter(outpath)
    logger.info('changing ID to url')

    try:
        df1['RSID'] = df1['RSID'].apply(lambda x: make_hyperlink(x, 'RSID'))
        df1['GENE_ID'] = df1['GENE_ID'].apply(lambda x: make_hyperlink(x, 'GENE'))
    except:
        logger.error('Cant parse ID Field')

    # temp column drop until applied in config
    df1.sort_index(inplace=True)
    df1.drop(columns=['DISTANCE', 'GENE_NAME', 'ERRORS / WARNINGS / INFO'], inplace=True)
    # reordering columns so ID and GENE ID are first
    cols_selected = ['GENE_ID', 'RSID', 'EFFECT', 'IMPACT', 'HGVS.P', 'HGVS.C', 'VARTYPE', '1000GP3_AF',
                     '1000GP3_AFR_AF', '1000GP3_AMR_AF', '1000GP3_EAS_AF', '1000GP3_EUR_AF', '1000GP3_SAS_AF',
                     'ESP6500_MAF_EA', 'ESP6500_MAF_AA', 'ESP6500_MAF_ALL', 'CLINVAR_CLNSIG', 'CLINVAR_CLNDSDB',
                     'CLINVAR_CLNDSDBID', 'CLINVAR_CLNDBN', 'CLINVAR_CLNREVSTAT', 'CLINVAR_CLNACC', 'ZIGOSITY',
                     'PolyPhen_Pred', 'PolyPhen_Score']
    # collist = [x for x in df1.columns if x not in cols_selected]

    df1 = df1[[x for x in cols_selected if x in df1.columns]]
    # df1 = df1.loc[:, df1.columns.isin(cols_selected)]
    # just for fleni, temp until i do it form cfg

    # removing qual column if duos or trios
    singlecols = ['QUAL']
    if any(["DUOS", "TRIOS"]) in df1.columns:
        df1.drop(columns=singlecols, inplace=True)
    workbook = output.book
    datasheet = workbook.add_worksheet('DATA')
    statsheet = workbook.add_worksheet('STATISTICS')
    output.sheets['DATA'] = datasheet
    format1 = workbook.add_format({'num_format': '###,###,###'})
    hyper = workbook.add_format({'hyperlink': True})
    datasheet.set_column('B:B', 18, format1)
    datasheet.set_column('D:D', 18, hyper)
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
        pass
    try:
        remove('./triostmp.png')
    except:
        pass

    try:
        remove('./duostmp.png')
    except:
        pass
