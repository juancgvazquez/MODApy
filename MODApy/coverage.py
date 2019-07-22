import logging
import multiprocessing as mp
import subprocess

import pandas as pd

from MODApy.cfg import cfg

logger = logging.getLogger(__name__)


def generate_coverage(file):
    cmd = "bedtools"
    args = " genomecov -ibam {} -bga".format(file)
    with open('{}_genomecov.bed'.format(file.rsplit('.', maxsplit=1)[0]), 'w') as output:
        p = subprocess.Popen(cmd + args, stdout=output, stderr=output, shell=True)
        p.wait()


def panel_intersect(file, panel_file):
    cmd = "bedtools"
    args = " intersect -a {} -b {}".format(file, panel_file)
    with open('{}_{}.bed'.format(file.rsplit('.', maxsplit=1)[0], panel_file.rsplit('.', maxsplit=1)[0]),
              'w') as output:
        subprocess.Popen(cmd + args, stdout=output, stderr=output, shell=True)


def annotate_genes(file, gene_file):
    cmd = "bedtools"
    args = " intersect -a {} -b {} -wb | awk -v OFS='\t' '{{print $1,$2,$3,$4,$8,$9}}'".format(file, gene_file)
    with open('{}_withgenes.cov'.format(file.rsplit('.', maxsplit=1)[0]), 'w') as output:
        p2 = subprocess.Popen(cmd + args, stdout=output, stderr=output, shell=True)
        p2.wait()


def create_coverage_reports(file):
    file = file.split('.')[0]
    csv = pd.read_csv('./' + file + '_genomecov_withgenes.cov', sep='\t',
                      names=['CHROM', 'START', 'END', 'DEPTH', 'GENE', 'STRAND'])
    csv['EXON'] = csv['GENE'].str.rsplit('-', n=1, expand=True)[1]
    csv['GENE'] = csv['GENE'].str.rsplit('-', n=1, expand=True)[0]
    csv.groupby(['GENE', 'CHROM']).DEPTH.describe().to_csv('./' + file + '_coverage_per_gene.csv')
    csv.groupby(['GENE', 'EXON', 'CHROM']).DEPTH.describe().to_csv('./' + file + '_coverage_per_gene.csv')


def main(files, bedfile, panelfile=None):
    pool = mp.Pool(processes=int(cfg['GENERAL']['cores']))
    try:
        logger.info('Generating coverage on {}'.format(files))
        r = pool.map_async(generate_coverage, files)
        r.wait()
        covfiles = ['{}_genomecov.bed'.format(file.rsplit('.', maxsplit=1)[0]) for file in files]
        intfiles = ['{}_{}.bed'.format(file.rsplit('.', maxsplit=1)[0], panelfile.rsplit('.', maxsplit=1)[0]) for file
                    in files]
        genfiles = ['{}_with_genes.cov'.format(file.rsplit('.', maxsplit=1)[0]) for file in files]
        if panelfile is not None:
            logger.info('Intersecting panel in files {}'.format(covfiles))
            panel = pool.starmap_async(panel_intersect, [(file, panelfile) for file in covfiles])
            panel.wait()
            logger.info('Annotating genes {}'.format(intfiles))
            ann = pool.starmap_async(annotate_genes, [(file, bedfile) for file in intfiles])
            ann.wait()
        else:
            logger.info('Annotating genes {}'.format(covfiles))
            ann = pool.starmap_async(annotate_genes, [(file, bedfile) for file in covfiles])
            ann.wait()
        logger.info('Creating Gene coverage report on {}'.format(genfiles))
        report = pool.map_async(create_coverage_reports, genfiles)
        report.wait()
        pool.close()
        pool.join()
    except Exception as e:
        logger.error('Error running coverage. Check logs for debugging.')
        logger.debug(e)
