import json
import logging
import os

import requests
from tqdm import tqdm

from MODApy import cfg

logger = logging.getLogger(__name__)
downlog = cfg.rootDir + '/logs/downloads.log'


def get_links(filename):
    link_list = list()
    if (filename.startswith('http://') or filename.startswith('https://') or filename.startswith('ftp://')):
        logger.debug('URL to download')
        download(filename)
    elif (filename.startswith('www')):
        logger.error('URLs must start with http:// or https:// or ftp://')
        exit(1)
    elif os.path.exists(filename):
        if filename.rsplit('.')[-1] == 'xlsx':
            logger.info('Parsing file to find URLs')
            from openpyxl import load_workbook
            wb = load_workbook(filename=filename)
            for ws in wb:
                for row in ws.rows:
                    for cell in row:
                        if cell.data_type == 'f':
                            link_list.append(cell.value.strip('=HYPERLINK("').split(',')[0])
            link_list = [x for x in link_list if '.tar' in x]
            logger.info('Parsing finished. Found %i links. These files will be now downloaded' % (len(link_list)))
            with open(downlog, 'r') as dlog:
                down_dict = json.load(dlog)
                down_dict.update({x: 'Pending' for x in link_list})
            with open(downlog, 'w') as dlog:
                json.dump(down_dict, dlog, indent=4)
            for link in link_list:
                download(link)
        elif filename.split('.')[-1] == 'xls':
            logger.info('Parsing file to find URLs')
            import re
            with open(filename, "r", encoding='ISO-8859-1') as fp:
                pru = fp.read()
            lista = re.findall(r'(https?://\S+)', pru)
            for x in lista:
                link_list.append(x.rsplit(sep='\x17')[0])
            link_list = [x for x in link_list if '.tar' in x]
            logger.info('Parsing finished. Found %i links. These files will be now downloaded' % (len(link_list)))
            with open(downlog, 'r') as dlog:
                down_dict = json.load(dlog)
                down_dict.update({x: 'Pending' for x in link_list})
            with open(downlog, 'w') as dlog:
                json.dump(down_dict, dlog, indent=4)
            for link in link_list:
                download(link)
        else:
            logger.error('File extension must be xlsx or xls.')
            exit(1)
    else:
        logger.error('File or URL provided is not valid')
        exit(1)


def download(url):
    """
    @param: url to download file
    """
    outputdir = cfg.patientPath + url.rsplit('/')[-1].split('.')[0] + '/'
    outputfile = url.rsplit('/')[-1]
    outpath = outputdir + outputfile
    tmpdir = cfg.rootDir + '/tmp/downloads/'
    tmppath = cfg.rootDir + '/tmp/downloads/' + outputfile
    os.makedirs(outputdir, exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)
    logger.info('Downloading %s to %s' % (url, outpath))
    try:
        response = requests.head(url)
    except:
        logger.error('Connection Failed')
        logger.debug('', exc_info=True)
        exit(1)
    else:
        if response.status_code == 200:
            if 'Content-Length' in response.headers.keys():
                file_size = int(response.headers["Content-Length"])
            else:
                logger.info('Could not get file size')
                file_size = 0
        else:
            logger.error('Connection Failed')
            logger.debug('', exc_info=True)
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
            logger.info('Moving file to destination folder')
            os.rename(tmppath, outpath)
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
            logger.error('', exc_info=True)
        with open(downlog, 'w') as dlog:
            json.dump(down_dict, dlog, indent=4)
    return file_size
