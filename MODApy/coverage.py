import multiprocessing as mp
import subprocess
import sys

import pandas as pd


def generate_coverage(file):
    cmd = "bedtools"
    args = " genomecov -ibam {} -bga".format(file)
    with open('{}_genomecov.bed'.format(file.split('.')[0]), 'w') as output:
        p = subprocess.Popen(cmd + args, stdout=output, stderr=output, shell=True)
        p.wait()


def panel_intersect(file, panel_file):
    cmd = "bedtools"
    args = " intersect -a {} -b {}".format(file, panel_file)
    with open('{}_{}.bed'.format(file.split('.')[0], panel_file.split('.')[0]), 'w') as output:
        subprocess.Popen(cmd + args, stdout=output, stderr=output, shell=True)


def annotate_genes(file, gene_file):
    cmd = "bedtools"
    args = " intersect -a {} -b {} -wb | awk -v OFS='\t' '{{print $1,$2,$3,$4,$8,$9}}'".format(file, gene_file)
    with open('{}_withgenes.cov'.format(file.split('.')[0]), 'w') as output:
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


if __name__ == '__main__':
    ##TODO Make argparser
    files_list = sys.argv[1:]
    b = sys.argv[1]
    pool = mp.Pool(processes=5)
    try:
        #		r= pool.map_async(generate_coverage,a)
        #		r.wait()
        pool.starmap(annotate_genes, [(file, b) for file in a])
        pool.close()
        pool.join()
    except Exception as e:
        print(str(e))
