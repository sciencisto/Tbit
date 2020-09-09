#__requires__ = "bcr-models==0.0.0"
#import re
#import sys
#import os

#import lyra packages
import bcr_models as bcr
import bcr_models.utils
import bcr_models.canonical_structures

#from pprint import pprint
#import re
#import csv


hmms = bcr.db.builtin_hmms()
template_db = bcr.db.BuiltinTemplateDatabase()
pdb_db = bcr.db.BuiltinPDBDatabase()
#csdb = bcr.db.BuiltinCsDatabase()

tmp_seq = "GDSVIQMQGQVTFSENDSLFINCTYSTTGYPTLFWYVQYSGEGPQLLLQVTTANNKGSSRGFEATYDKGTTSFHLQKTSVQEIDSAVYYCAISDLSGGSNAKLAFGKGTKLSVK"
ig_chain =bcr.IgChain(tmp_seq, template_db=template_db, pdb_db=pdb_db)
#Score

ig_chain.hmmsearch(*hmms)
if (ig_chain.E_value<1e-15 ):
    print(ig_chain)
    print(ig_chain.aligned_seq)
