import logging
from os import path, remove

import matplotlib
import pandas as pd

from MODApy.cfg import cfg
from MODApy.vcfmgr import ParsedVCF

matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib_venn as venn

logger = logging.getLogger(__name__)

'''
Helper function to get statistics out of vcfs
'''


def getstats(self, type=0):
    stats = {}
    if (type == 0):
        if all(col in self.columns for col in ['ZIGOSITY', 'VARTYPE', 'IMPACT', 'EFFECT']):
            vcfstats = self.groupby(['CHROM', 'ZIGOSITY', 'VARTYPE', 'IMPACT',
                                     'EFFECT']).size().to_frame(name='count')
            vcfstats.name = 'stats'
            stats['df'] = []
            stats['df'].append(vcfstats)
            plt.pie(vcfstats.groupby('CHROM').count(), labels=vcfstats.groupby('CHROM').size().index.values)
            my_circle = plt.Circle((0, 0), 0.7, color='white')
            chromVars = plt.gcf()
            chromVars.gca().add_artist(my_circle)
            stats['graphs'] = []
            stats['graphs'].append(chromVars)
    elif (type == 1):
        if 'TRIOS' in self.columns:
            trios = self.groupby('TRIOS', sort=False).size()
            names = self.name.split(':')
            names.append(self.name)
            names.append(names[0] + ':' + names[1])
            names.append(names[1] + ':' + names[2])
            names.append(names[0] + ':' + names[2])
            trios = trios.reindex(names).fillna(0)
            trios = trios.astype(int)
            A, B, C = self.name.split(':')
            triosgraph = plt.figure()
            venn.venn3(trios, set_labels=[A, B, C], set_colors=['b', 'r', 'g'])
            triosgraph.savefig('./triostmp.png', dpi=triosgraph.dpi)
            triosgraph.clf()
        elif 'DUOS' in self.columns:
            duos = self.groupby('DUOS', sort=False).size()
            names = self.name.split(':')
            names.append(self.name)
            duos = duos.reindex(names).fillna(0)
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

    logger.error("%s couldn't be found. Please check if file exists and that it's extension is %s" % (filePath,
                                                                                                      extension))
    exit(1)


def aminoChange(value: str):
    value = value.replace('p.', '')
    if value[:3] != value[-3:]:
        return 'CHANGE'


def df_to_excel(df1: ParsedVCF, outpath):
    #    output = pd.ExcelWriter(outpath, engine='xlsxwriter', options={'strings_to_urls': False})
    output = pd.ExcelWriter(outpath)
    # temp column drop until applied in config
    # df1.reset_index(inplace=True)

    # reordering columns according to cfg
    cols_selected = cfg["OUTPUT"]["columnsorder"].replace(',', ' ').split()

    df1['AMINOCHANGE'] = df1['HGVS.P'].apply(aminoChange)

    df1 = df1[[x for x in cols_selected if x in df1.columns]].copy()
    df1 = df1.sort_values(by=cols_selected[0])

    workbook = output.book
    datasheet = workbook.add_worksheet('DATA')
    statsheet = workbook.add_worksheet('STATISTICS')

    output.sheets['DATA'] = datasheet
    if ('CLINVAR_CLNSIG' in df1.columns):
        translation = {'255': 'other', '0': 'Uncertain significance', '1': 'not provided', '2': 'Benign',
                       '3': 'Likely Benign', '4': 'Likely pathogenic', '5': 'Pathogenic', '6': 'drug response',
                       '7': 'histocompatibility'}
        for k, v in translation.items():
            df1['CLINVAR_CLNSIG'] = df1['CLINVAR_CLNSIG'].str.replace(k, v)

    formatnum = workbook.add_format({'num_format': '#,#####0.00000'})
    for i, col in enumerate(df1.columns):
        # collen = df1[col].astype(str).map(len).max()
        # collen = max(collen, len(col) + 4)
        datasheet.set_column(i, i, 15, formatnum)

    formatpos = workbook.add_format({'num_format': '###,###,###'})
    datasheet.set_column(cols_selected.index('POS') - 2, cols_selected.index('POS') - 2, 15, formatpos)
    datasheet.set_column(cols_selected.index('RSID'), cols_selected.index('RSID'), 15)
    # Light red fill with dark red text.
    highformat = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006', 'bold': True})
    # Light yellow fill with dark yellow text.
    modformat = workbook.add_format({'bg_color': '#FFFF99', 'font_color': '#9C6500', 'bold': True})
    # Light orange fill with dark orange text.
    moderformat = workbook.add_format({'bg_color': '#FFCC99', 'font_color': '#FF6600', 'bold': True})
    # Green fill with dark green text.
    lowformat = workbook.add_format({'bg_color': '#C6EFCE', 'font_color': '#006100', 'bold': True})
    datasheet.conditional_format(0, cols_selected.index('IMPACT'),
                                 len(df1), cols_selected.index('IMPACT'),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'HIGH', 'format': highformat, })
    datasheet.conditional_format(0, cols_selected.index('IMPACT'),
                                 len(df1), cols_selected.index('IMPACT'),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'MODIFIER', 'format': modformat})
    datasheet.conditional_format(0, cols_selected.index('IMPACT'),
                                 len(df1), cols_selected.index('IMPACT'),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'MODERATE', 'format': moderformat})
    datasheet.conditional_format(0, cols_selected.index('IMPACT'),
                                 len(df1), cols_selected.index('IMPACT'),
                                 {'type': 'text', 'criteria': 'containing', 'value': 'LOW', 'format': lowformat})
    df1.to_excel(output, sheet_name='DATA', merge_cells=False, index=False)

    if (df1.reset_index().index.max() < 32150):
        logger.info('changing ID to url')
        try:
            colid = cols_selected.index('RSID')
            colgen = cols_selected.index('GENE_NAME')
            row = 2
            for x in zip(df1['RSID'], df1['GENE_NAME']):
                if type(x[0]) == str:
                    urlrs = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=%s"
                    rsvalue = (x[0].replace(';', ',').split(','))[0]
                    datasheet.write_url('%s%i' % (chr(colid + 65), (row)),
                                        urlrs % rsvalue, string=rsvalue)
                if type(x[1]) == str:
                    urlgen = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s"
                    datasheet.write_url('%s%i' % (chr(colgen + 65), (row)),
                                        urlgen % x[1], string=x[1])
                row += 1
        except Exception as e:
            logger.error(e, exc_info=True)
    datasheet.autofilter(0, 0, len(df1), len(cols_selected))
    stats = getstats(df1)
    output.sheets['STATISTICS'] = statsheet
    try:
        for i in range(len(stats['df'])):
            stats['df'][i].to_excel(output, sheet_name='STATISTICS')
    except Exception as e:
        logger.error('Could not print statistics. Error was {}'.format(e), exc_info=True)
    try:
        for i in range(len(stats['graphs'])):
            stats['graphs'][i].savefig('./tempgraph.png')
            statsheet.insert_image('H2', './tempgraph.png')
    except Exception as e:
        logger.error('Could not print stats graphs. Error was {}'.format(e), exc_info=True)
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
