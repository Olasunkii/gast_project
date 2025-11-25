import sys
from sra_extractor import SRAExtractor
from Bio import Entrez

email=sys.argv[1]
organism=sys.argv[2]
retmax=int(sys.argv[3])

Entrez.email=email
ex=SRAExtractor()
# Fetch metadata and AST table
df=ex.resistant_paired_metadata(organism,retmax=retmax)
fname=f"SraRunInfo_{organism.replace(' ','_')}.csv"# Save metadata)
ex.save_metadata(df,fname)
