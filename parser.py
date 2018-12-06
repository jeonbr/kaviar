import os
import csv
from biothings.utils.dataload import unlist
from biothings.utils.dataload import value_convert_to_number
from biothings.utils.dataload import merge_duplicate_rows, dict_sweep
from utils.hgvs import get_hgvs_from_vcf
from itertools import groupby, chain
import json
from tempfile import mkstemp
import subprocess
import vcf
import tarfile
import gzip, shutil

VALID_COLUMN_NO = 8

'''this parser is for Kaviar version 160204-Public(All variants, annotated with data sources) downloaded from
http://db.systemsbiology.net/kaviar/Kaviar.downloads.html'''


# convert one snp to json
# convert one snp to json
def _map_line_to_json(item):
    chrom = item.CHROM
    chromStart = item.POS
    ref = item.REF
    info = item.INFO

    try:
        af = info['AF']
    except:
        af = None
    try:
        ac = info['AC']
    except:
        ac = None
    try:
        an = info['AN']
    except:
        ac = None
    try:
        ds = info['DS']
    except:
        ds = None

    # convert vcf object to string
    item.ALT = [str(alt) for alt in item.ALT]

    # if multiallelic, put all variants as a list in multi-allelic field
    hgvs_list = None
    if len(item.ALT) > 1:
        hgvs_list = []
        for alt in item.ALT:
            try:
                hgvs_list.append(get_hgvs_from_vcf(chrom, chromStart, ref, alt, mutant_type=False))
            except:
                hgvs_list.append(alt)

        assert len(item.ALT) == len(info['AC']), "Expecting length of item.ALT= length of info.AC, but not for %s" % (item)
        assert len(item.ALT) == len(info['AF']), "Expecting length of item.ALT= length of info.AF, but not for %s" % (item)
        if ds:
            if len(item.ALT) != len(info['DS']):
                ds_str = ",".join(info['DS'])
                ds_str = ds_str.replace("NA7022,18", "NA7022_18")
                ds_list = ds_str.split(",")
                info['DS'] = [d.replace("NA7022_18", "NA7022,18") for d in ds_list]
                assert len(item.ALT) ==len(info['DS']), "info.DS mismatch, %s: %s\n## DS: %s" % (item, info['DS'])

    for i, alt in enumerate(item.ALT):
        try:
            (HGVS, var_type) = get_hgvs_from_vcf(chrom, chromStart, ref, alt, mutant_type=True)
        except:
            continue

        if HGVS is None:
            return

        # load as json data
        one_snp_json = {
            "_id": HGVS,
            "kaviar": {
                "multi-allelic": hgvs_list,
                "ref": ref,
                "alt": alt,
                "af": info['AF'][i],
                "ac": info['AC'][i],
                "an": an,
                "ds": info['DS'][i].split("|") if ds else None,
            }
        }

        yield value_convert_to_number(one_snp_json)



# open file, parse, pass to json mapper
def load_data(data_folder):
    tar = tarfile.open(os.path.join(data_folder, "Kaviar-160204-Public-hg19.vcf.tar"))
    tar.extractall(data_folder)
    tar.close()
    with gzip.open(os.path.join(data_folder, "Kaviar-160204-Public", "vcfs", "Kaviar-160204-Public-hg19.vcf.gz"), 'r' ) as f_in:
        with open(os.path.join(data_folder, "Kaviar-160204-Public-hg19.vcf"), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)    
                                         
    input_fn = os.path.join(data_folder,"Kaviar-160204-Public-hg19.vcf")
    vcf_reader = vcf.Reader(open(input_fn, 'r'), strict_whitespace=True)
    json_rows = map(_map_line_to_json, vcf_reader)
    json_rows = chain.from_iterable(json_rows)

    fd_before_sort, temp_path_before_sort = mkstemp(dir=data_folder)
    fd_after_sort, temp_path_after_sort = mkstemp(dir=data_folder)    

    with open(temp_path_before_sort, "w") as f:
        dbwriter = csv.writer(f)
        for doc in json_rows:
            dbwriter.writerow([doc['_id'], json.dumps(doc)]) 

    popen_str = 'sort {} > {}'.format(temp_path_before_sort, temp_path_after_sort)
    p = subprocess.Popen(popen_str, shell=True)
    os.waitpid(p.pid, 0)
    os.close(fd_before_sort)
    
    json_rows = csv.reader(open(temp_path_after_sort))
    json_rows = (json.loads(row[1]) for row in json_rows)
    row_groups = (it for (key, it) in groupby(json_rows, lambda row: row["_id"]))
    json_rows = (merge_duplicate_rows(rg, "kaviar") for rg in row_groups)
    return (unlist(dict_sweep(row, vals=[None, ])) for row in json_rows)


