import logging
import multiprocessing as mp
import subprocess
import traceback

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
    args = " intersect -a {} -b {} -wb | awk -v OFS='\t' '{{print $1,$2,$3,$4,$8}}'".format(file, gene_file)
    with open('{}_with_genes.cov'.format(file.rsplit('.', maxsplit=1)[0]), 'w') as output:
        p2 = subprocess.Popen(cmd + args, stdout=output, stderr=output, shell=True)
        p2.wait()


def create_coverage_reports(file):
    # file = file.rsplit('.', maxsplit=1)[0]
    logger.debug(file)
    csv = pd.read_csv(file, sep='\t', names=['CHROM', 'START', 'END', 'DEPTH', 'EXON'])
    csv['GENE'] = csv['EXON'].str.split('_', n=1).str[0]
    csv.groupby(['GENE', 'CHROM']).DEPTH.describe().drop(columns='count').to_csv(file.rsplit('.',maxsplit=1)[0] + '_coverage_per_gene.csv')
    csv.groupby(['EXON', 'CHROM']).DEPTH.describe().drop(columns='count').to_csv(file.rsplit('.',maxsplit=1)[0] + '_coverage_per_exon.csv')


def main(files, bedfile, panelfile=None):
    pool = mp.Pool(processes=int(cfg['GENERAL']['cores']))
    try:
        logger.info('Generating coverage on {}'.format(files))
        r = pool.map_async(generate_coverage, files)
        r.wait()
        covfiles = ['{}_genomecov.bed'.format(file.rsplit('.', maxsplit=1)[0]) for file in files]
        genfiles = ['{}_genomecov_with_genes.cov'.format(file.rsplit('.', maxsplit=1)[0]) for file in files]
        if panelfile is not None:
            logger.info('Intersecting panel in files {}'.format(covfiles))
            panel = pool.starmap_async(panel_intersect, [(file, panelfile) for file in covfiles])
            panel.wait()
            intfiles = ['{}_{}.bed'.format(file.rsplit('.', maxsplit=1)[0], panelfile.rsplit('.', maxsplit=1)[0]) for
                        file
                        in files]
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
        logger.debug(traceback.print_exc())
