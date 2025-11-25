import sys
from sra_extractor import SRAExtractor
from Bio import Entrez

email=sys.argv[1]
organism=sys.argv[2]
retmax=int(sys.argv[3])

Entrez.email=email
ex=SRAExtractor()
# Fetch metadata
df=ex.fetch_runinfo_paired(organism,retmax=retmax)
# Save enriched metadata (so we keep size_human, reads, etc.)
fname=f"SraRunInfo_{organism.replace(' ','_')}.csv"
ex.save_metadata(df,fname)
