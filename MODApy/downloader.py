import hashlib
import json
import logging
import os
import re
import tarfile

import openpyxl
import pandas as pd
import requests
import xlrd
from tqdm import tqdm

from MODApy import cfg

logger = logging.getLogger(__name__)
downlog = cfg.rootDir + '/logs/downloads.log'


def get_links(filename):
    link_list = list()
    # If filename is an url
    if (filename.startswith('http://') or filename.startswith('https://') or filename.startswith('ftp://')):
        logger.debug('URL to download')
        download(filename)
    elif (filename.startswith('www')):
        logger.error('URLs must start with http:// or https:// or ftp://')
        exit(1)
    # If filename is an excel file
    elif os.path.exists(filename):
        if filename.rsplit('.')[-1] == 'xlsx':
            logger.info('Parsing file to find URLs')
            wb = openpyxl.load_workbook(filename=filename)
            for ws in wb:
                for row in ws.rows:
                    for cell in row:
                        if cell.data_type == 'f':
                            link_list.append(cell.value.strip('=HYPERLINK("').replace('"', '').split(',')[0])
            pdxl = pd.read_excel(filename, sheet_name='Download_Address')
            md5 = list(pdxl.iloc[pdxl.index[pdxl['download address'] == 'md5sum'].item() + 1:]['download address'])
            link_list = [x for x in link_list if '.tar' in x]
            linksdict = dict(zip(link_list, md5))
            logger.info('Parsing finished. Found %i links. These files will be now downloaded' % (len(link_list)))
            with open(downlog, 'r') as dlog:
                down_dict = json.load(dlog)
                down_dict.update({x: 'Pending' for x in linksdict})
            with open(downlog, 'w') as dlog:
                json.dump(down_dict, dlog, indent=4)
            for link in linksdict:
                download(link, linksdict[link])
        elif filename.split('.')[-1] == 'xls':
            wb = xlrd.open_workbook(filename)
            ws = wb.sheet_by_name('Download_Address')
            for row in range(ws.nrows):
                for col in range(ws.ncols):
                    if ws.cell(row, col).value == 'md5sum':
                        md5list = ws.col_values(col, row + 1)
            logger.info('Parsing file to find URLs')
            with open(filename, "r", encoding='ISO-8859-1') as fp:
                pru = fp.read()
            lista = re.findall(r'(https?://\S+)', pru)
            for x in lista:
                link_list.append(x.rsplit(sep='\x17')[0])
            link_list = [x for x in link_list if '.tar' in x]
            linksdict = dict(zip(link_list, md5list))
            logger.info('Parsing finished. Found %i links. These files will be now downloaded' % (len(link_list)))
            with open(downlog, 'r') as dlog:
                down_dict = json.load(dlog)
                down_dict.update({x: 'Pending' for x in linksdict})
            with open(downlog, 'w') as dlog:
                json.dump(down_dict, dlog, indent=4)
            for link in linksdict:
                download(link, linksdict[link])
        else:
            logger.error('File extension must be xlsx or xls.')
            exit(1)
    else:
        logger.error('File or URL provided is not valid')
        exit(1)


def download(url, md5=None):
    """
    @param: url to download file
    """
    outputdir = cfg.patientPath + url.rsplit('/')[-1].split('.')[0] + '/'
    outputfile = url.rsplit('/')[-1]
    outpath = outputdir + outputfile
    tmpdir = cfg.tmpDir
    tmppath = cfg.tmpDir + outputfile
    os.makedirs(outputdir, exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)
    logger.info('Downloading %s to %s' % (url, outpath))
    try:
        response = requests.head(url)
    except:
        logger.info('Connection Failed')
        logger.debug('Downloader Connection Failed', exc_info=True)
        exit(1)
    else:
        if response.status_code == 200:
            if 'Content-Length' in response.headers.keys():
                file_size = int(response.headers["Content-Length"])
            else:
                logger.info('Could not get file size')
                file_size = 0
        else:
            logger.info('Connection Failed')
            logger.debug('Downloader Connection Failed', exc_info=True)
            exit(1)
    if os.path.exists(tmppath):
        first_byte = os.path.getsize(tmppath)
    else:
        first_byte = 0
    header = {"Range": "bytes=%s-%s" % (first_byte, file_size)}
    with open(downlog, 'r') as dlog:
        down_dict = json.load(dlog)
        down_dict.update({url: 'Downloading'})
    with open(downlog, 'w') as dlog:
        json.dump(down_dict, dlog, indent=4)
    pbar = tqdm(
        total=file_size, initial=first_byte,
        unit='B', unit_scale=True, desc=url.split('/')[-1])
    try:
        req = requests.get(url, headers=header, stream=True)
        with(open(tmppath, 'ab')) as f:
            for chunk in req.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
                    pbar.update(1024)
        pbar.close()
        with open(downlog, 'r') as dlog:
            down_dict = json.load(dlog)
            down_dict.update({url: 'Download Complete'})
            logger.info('Download Complete for %s' % outputfile)
            if md5 is not None:
                logger.info('Checking md5')
                file_md5 = hashlib.md5()
                with open(tmppath, 'rb') as f:
                    while True:
                        fileBuffer = f.read(16 * 1024 * 1024)
                        if not fileBuffer:
                            break
                        file_md5.update(fileBuffer)
                    file_md5 = file_md5.hexdigest()
                if md5 == file_md5:
                    logger.info('MD5 Verified')
                    logger.info('Moving file to destination folder')
                    os.rename(tmppath, outpath)
                    try:
                        logger.info('Extracting tar file')
                        tf = tarfile.open(outpath)
                        tf.extractall(outputdir)
                        logger.info('Extraction finished')
                    except:
                        logger.error('Extraction failed')
                        logger.debug('', exc_info=True)
                        exit(1)
                else:
                    print(file_md5, md5)
                    logger.info('MD5 Verification Failed')
                    logger.error('Error with the download')
                    exit(1)
            try:
                os.removedirs(tmpdir)
            except:
                logger.debug(tmpdir + ' is not empty.')
        with open(downlog, 'w') as dlog:
            json.dump(down_dict, dlog, indent=4)
    except:
        with open(downlog, 'r') as dlog:
            down_dict = json.load(dlog)
            down_dict.update({url: 'Download Error'})
            logger.error('Error with download')
            logger.debug('', exc_info=True)
        with open(downlog, 'w') as dlog:
            json.dump(down_dict, dlog, indent=4)
    return file_size
