"""
This script should be executed once, before any .vcf annotations.
Runtime exceeds 15 minutes because it fetches all unique accession ids
from clinVar.

It will save all "pathogenic" and/or "conflicting" SNVs and Indels as a
dictionary object in a pickle file in the current working directory.

After the initial execution of clinVar_parser.py it's uneccessary to run it
again, until clinVar release a new update, or whenever you'd like.
"""
import os, re, time
from xml.dom import minidom

try:
    import json
except ImportError:
    import simplejson as json

try: #Python2.x
    from urllib2 import urlopen
    from urllib import urlencode
except ImportError: #Python3.x
    from urllib.request import urlopen
    from urllib.parse import urlencode

try:
    import cPickle as pickle
except ImportError:
    import cPickle

def esearch(chrom):
    """
    Returns only pathogenic or conflicting pathogenic variants from
    the NCBI Entrez EUtils API
    """   
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    # includes conflicting variants
    esearch_url = base_url + "esearch.fcgi?db=clinvar&term=%s[chr]+AND+clinsig_pathogenic[prop]&retmode=json&retmax=50000"
    url =  esearch_url % (chrom)
    response = urlopen(url)
    search_data = json.loads(response.read())
    
    return search_data

def efetch(chr_uids):
    """
    Returns a dictionary of each pathogenic clinVar accession with HGVS:
        
    {NP_001270035.1: {'224774': {'GRCh38': 'g.20043503C>T', 'AAconseq': 'p.Arg140Ter', 'GRCh37': 'g.20043503C>T'}})
    
    Currently returns 33,268 clinVar accessions of the ~47,000 accessions.
    The remaining 13,760 are Indels and copy number variants where there is no HGVS protein change.
    """
    chr_uids = [str(uid) for uid in chr_uids['esearchresult']['idlist']]
    uidstr = ','.join(chr_uids) # stringify the list of chr_uids
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    efetch_url = base_url + "efetch.fcgi"
    params = urlencode({'db': 'clinvar',
              'rettype': 'variation',
              'id': uidstr})
    response = urlopen(efetch_url, data=params) #POST per efetch guidelines
    search_data = minidom.parse(response)
    
    d = {}
    for report in search_data.getElementsByTagName('VariationReport'):
        uid = report.attributes['VariationID'].value
        # each uid can have multiple of each.
        RefSeq = []
        AAconseq = []
        GRCh38 = []
        GRCh37 = []
        molConseq = []
        
        try:
            for elem in report.getElementsByTagName('HGVS'):
                if "HGVS, protein, RefSeq" in elem.attributes['Type'].value:
                    RefSeq.append(elem.attributes['AccessionVersion'].value)
                    AAconseq.append(elem.attributes['Change'].value)
                
                if "HGVS, genomic, top level" == elem.attributes['Type'].value:
                    GRCh38.append(elem.attributes['Change'].value)
                
                if "HGVS, genomic, top level, previous" == elem.attributes['Type'].value:
                    GRCh37.append(elem.attributes['Change'].value)
        
        except KeyError as e:
            # some (<10) uids have HGVS without 'Change' attribute
            continue
        
        try:
            for elem in report.getElementsByTagName('MolecularConsequence'):
                    molConseq.append(elem.attributes['Function'].value)
        except KeyError as e:
            continue
             
        # Create nested dictionary -- keys are RefSeq transcript id.
        for rseq,aa in zip(RefSeq,AAconseq):
            try:
                d[rseq][uid] = {"AAconseq": aa,
                                "GRCh38": list(set(GRCh38)),
                                "GRCh37": list(set(GRCh37)),
                                "molConseq": list(set(molConseq))
                               }
            except KeyError as e: #no entry, create new
                d[rseq] = {}
                d[rseq][uid] = {"AAconseq": aa,
                                "GRCh38": list(set(GRCh38)),
                                "GRCh37": list(set(GRCh37)),
                                "molConseq": list(set(molConseq))
                               }
    
    return d

def merge_dicts(*dict_args):
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

if __name__ == '__main__':
    chromstr = ['X','Y','MT']
    chromint = [x for x in xrange(1,23)]
    chroms = chromint + chromstr
    
    # purposely NOT using multiprocessing here, don't want to flood NCBI
    # additionaly - time.sleep(1) - NCBI recommends 3 requests/sec
    for chrom in chroms:
        chr_uids = esearch(chrom)
        cvObj = efetch(chr_uids)
        pkl_out = open("clinVar_chr%s.p" % chrom, 'wb')
        pickle.dump(cvObj, pkl_out)
        pkl_out.close()
        print ("Finished chromosome: %s" % chrom)
        time.sleep(30)
    
    # Iterating again because holding all the chroms in memory and trying to
    # process the last few chroms requires 5 GB+ RAM, so decided to save
    # each chrom separately as a pickle.
    fetchresults=[]
    for chrom in chroms:
        pkl_in = open('clinVar_chr%s.p' % chrom, 'rb')
        clinVar_obj = pickle.load(pkl_in)
        fetchresults.append(clinVar_obj)
        pkl_in.close()
        os.remove('clinvar_chr%s.p' % chrom)
    
    
    clinVar_out = merge_dicts(*fetchresults)
    output = open('clinVar_obj.p', 'wb')
    pickle.dump(clinVar_out, output)
    output.close()