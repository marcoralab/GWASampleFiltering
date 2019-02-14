import re
import glob
import shutil
import os
from snakemake.io import *


def fix_path(path):
    path = os.path.abspath(path)
    assert os.path.isdir(path), "Path not found: {}".format(path)
    return path


def build_samp(in_path, extension='bim', flat=True, excl=None, incl=None):
    p = os.path.abspath(in_path)
    if flat:
        glterm = p + "/*." + extension
        files = glob.iglob(glterm)
        rm = len(extension) + 1
        samples = [os.path.basename(x)[:-rm] for x in files]
    else:
        p = [directory for directory, y, files in os.walk(p)
             if any('.' + extension in f for f in files)]
        samples = [os.path.basename(x) for x in p]
    if excl:
        samples = [x for x in samples if excl not in x]
    if incl:
        samples = [x for x in samples if incl in x]
    return samples


def parser(config):
    BPLINK = ["bed", "bim", "fam"]
    DATAIN = fix_path(config['DataIn'])
    DATAOUT = fix_path(config['DataOut'])
    sample_conf = config['sample']
    FAMILY = 'T' if config['family'] else 'F'

    try:
        SEXIN = fix_path(config['SexIn'])
        sex_dir = True
    except:
        SEXIN = config['SexIn']
        sex_dir = False

    flat = len(build_samp(DATAIN, 'bim')) > 0

    if config['SexIn']:
        SAMPLE = build_samp(DATAIN, 'bim', flat, config['SexIn'])
        test_sex = build_samp(DATAIN, 'bim', flat, incl=config['SexIn'])
        sex_stem = len(test_sex) > 0
    else:
        SAMPLE = build_samp(DATAIN, 'bim', flat)

    if sample_conf:
        if type(sample_conf) == str:
            sample_conf = [sample_conf]
        diff = set(sample_conf) - set(SAMPLE)
        errstr = "The following samples cannot be found: "
        assert len(diff) == 0, errstr.format(diff)
        SAMPLE = list(set(sample_conf) & set(SAMPLE))

    start = {}

    if flat:
        start['files'] = expand("{DataIn}/{{sample}}.{ext}",
                                ext=BPLINK, DataIn=DATAIN)
        start['stem'] = expand("{DataIn}/{{sample}}", DataIn=DATAIN)[0]
        if config['SexIn']:
            assert sex_dir ^ sex_stem, "SexIn must be EITHER a stem OR dir"
            if sex_stem:
                start['sex'] = expand("{DataIn}/{{sample}}{SexIn}.{ext}",
                                      ext=BPLINK, DataIn=DATAIN, SexIn=SEXIN)
                start['sex_stem'] = expand("{DataIn}/{{sample}}{SexIn}",
                                           DataIn=DATAIN, SexIn=SEXIN)[0]
            elif sex_dir:
                start['sex'] = expand("{SexIn}/{{sample}}.{ext}",
                                      ext=BPLINK, SexIn=SEXIN)
                start['sex_stem'] = expand("{SexIn}/{{sample}}",
                                           SexIn=SEXIN)[0]
        else:
            start['sex'] = start['files']
            start['sex_stem'] = start['stem']
    else:
        start['files'] = expand("{DataIn}/{{sample}}/{{sample}}.{ext}",
                                ext=BPLINK, DataIn=DATAIN)
        start['stem'] = expand("{DataIn}/{{sample}}/{{sample}}",
                               DataIn=DATAIN)[0]
        if config['SexIn']:
            assert sex_dir ^ sex_stem, "SexIn must be EITHER a stem OR dir"
            if sex_dir:
                start['sex'] = expand("{SexIn}/{{sample}}/{{sample}}.{ext}",
                                      ext=BPLINK, SexIn=SEXIN)
                start['sex_stem'] = expand("{SexIn}/{{sample}}/{{sample}}",
                                           SexIn=SEXIN)[0]
            elif sex_stem:
                start['sex'] = expand(
                    "{DataIn}/{{sample}}/{{sample}}{SexIn}.{ext}",
                    ext=BPLINK, DataIn=DATAIN, SexIn=SEXIN)
                start['sex_stem'] = expand(
                    "{DataIn}/{{sample}}/{{sample}}{SexIn}",
                    DataIn=DATAIN, SexIn=SEXIN)[0]
        else:
            start['sex'] = start['files']
            start['sex_stem'] = start['stem']

    return start, FAMILY, SAMPLE, DATAOUT
