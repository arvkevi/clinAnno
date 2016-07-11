"""
The goal of this script is to annotate a .vcf or .vcf.gz file
when a variant occurs at the same amino acid position as a known
pathogenic variant in the clinVar database;
as defined in the guidelines developed by the ACMG:

Richards, Sue, et al.
"Standards and guidelines for the interpretation of sequence variants:
a joint consensus recommendation of the American College of
Medical Genetics and Genomics and the Association for Molecular Pathology."
Genetics in Medicine (2015).

Scope is currently limited to variants designated PM5 and/or PS1.
"""
import gzip
import re
import argparse

try:
    import cPickle as pickle
except ImportError:
    import cPickle

def process_variant(cvObj, record):
    """
    Function to handle one record of the .vcf
    Addtional functionality could be built here.
    """ 
    if record.startswith('#'):
        return record
    
    record = record.rstrip().split('\t')
    pm5_ps1_anno = pm5_ps1(record, cvObj) 
    if pm5_ps1_anno is None: # No annotations, write original vcf record
        return '\t'.join(record) + "\n"
    else:
        # if not None, output is a string:
        # 'PM5=5678;PS1=1234;5678' 
        return pm5_ps1_anno
    
def pm5_ps1(record, cvObj):
    """  
    The "PM5" desigsignation is applied when a variant is in the same
    AA residue as another known pathogenic variant, but the alternate
    Amino Acid differs from the previously reported variant.
    Example:
    KRT14 AA residue 119
    
    
    "The "PS1" designation is applied when a variant causes
    the same Amino Acid change regardless of nucleotide change in a residue
    that has previously been clasified as pathogenic.
    Example:
    ARID1B  AA residue 1718
    """
    
    record_string = '\t'.join(record)
    
    if "p." in record_string:
        
        # VEP uses this for synonymous
        if "p.%3D" in record_string:
            return
        
        # matches the first AA change (p.Lys119Arg)
        prot= re.match(r".*(NP_\d+.\d+):p.([a-zA-Z_?-]+)(\d+)([\w+_?-]{8,10}|[a-zA-Z_?-]{3}|=|\?)([a-zA-Z]+)?",
                       record_string)
        
        try:
            refseqNP = prot.group(1)
            AAref = prot.group(2)
            AApos = prot.group(3)
            AAalt = prot.group(4)
        
        ## TO DO: these are start_lost variants.
        except AttributeError:
            return
        
        #check for "fs" variant 
        try:
            AAfs = prot.group(5)
        except AttributeError:
            AAfs = None
        
        chrom = record[0]
        pos = record[1]
    
    else: #no AA change, skips record
        return
    
    # Dictionary that will be returned for annotating the vcf
    annoPath = {"PS1": [],
                "PM5": []
                    }
    
    # is the transcript in our clinVar object?
    try:
        clinVar_obj[refseqNP]
    except: #no
        return
    
    for uid, variant in clinVar_obj[refseqNP].items():
        try:
            prot= re.match(r".*p.([a-zA-Z_?-]+)(\d+)([\w+_?-]{8,10}|[a-zA-Z_?-]{3}|=|\?)([a-zA-Z]+)?",
                       str(variant["AAconseq"]))
            cvAAref = prot.group(1)
            cvAApos = prot.group(2)
            cvAAalt = prot.group(3)
        except:  # Gln44= was confusing, prob OK to remove try-except here.
            print "skipping uid: %s" % uid
            continue
        
        # match ref AA and AA position
        if AAref == cvAAref and AApos == cvAApos:
            if AAalt == cvAAalt: #PS1
                annoPath["PS1"].append(uid)
            
            else: #PM5
                if 'missense' in '\t'.join(clinVar_obj[refseqNP][uid]['molConseq']):
                    annoPath["PM5"].append(uid)
                    
    # don't write annotation on empty list
    if len(annoPath["PM5"]) == 0:
        del annoPath["PM5"]
    if len(annoPath["PS1"]) ==  0:
        del annoPath["PS1"]
    
    if annoPath == {}:
        return
    
    # translate the dictionary to a string
    pm5_ps1_str = ';'.join(["%s=" % (k) +
                            ';'.join(str(e) for e in v)
                            for k,v in annoPath.items()]) + ';'
    
    #prepend the annotation in the INFO column of the record
    record[7] = pm5_ps1_str + record[7]
    
    #return the string version of the record.
    return "\t".join(record) + "\n"

if __name__ == '__main__':
    parser = argparse.ArgumentParser()   
    parser.add_argument("--vcf_in", help="full path to input .vcf file")
    parser.add_argument("--vcf_out", help="full path to output .vcf file")
    parser.add_argument("--nproc", default=1, help="number of processes")
    
    args = parser.parse_args()

    f = args.vcf_in
    out = args.vcf_out
    
    try:
        clinvar_pkl = open('clinVar_obj.p', 'rb')
        clinVar_obj = pickle.load(clinvar_pkl)
        clinvar_pkl.close()
    except IOError as e:
        print "Run clinAnno_parser.py and make sure the .p file is in cwd"
        print e
    
    
    if f.endswith(".gz"):
        with gzip.open(f, 'rb') as v:
            with open (out, "w") as annoVCF:
                for line in v:
                    anno_record = process_variant(clinVar_obj, line)
                    annoVCF.write(anno_record)
    else:
        with open(f, 'rb') as v:
            with open (out, "w") as annoVCF:
                for line in v:
                    anno_record = process_variant(clinVar_obj, line)
                    annoVCF.write(anno_record)

