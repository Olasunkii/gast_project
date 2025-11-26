import sys
from sra_extractor import SRAExtractor
from Bio import Entrez

email=sys.argv[1]
organism=sys.argv[2]
retmax=int(sys.argv[3])

Entrez.email=email
ex=SRAExtractor()
# Fetch metadata and save corresponding AST table (to paht:data/ast)
metadata_df=ex.collect_resistant_metadata(organism,retmax=retmax)
fname=f"SraRunInfo_{organism.replace(' ','_')}.csv"# Save metadata)
ex.save_metadata(metadata_df,fname)#save metadata to csv file: data/metadata
