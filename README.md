# gast_project


### Prerequisites

* Install required packages for automated data retrieval.
* Download the associated databases (Checkm2 and bakta).


---
# Still work in progress

**Needed reference databases installation**
These needs to be installed before workflow running

install bakta_db:
wget https://zenodo.org/record/14916843/files/db-light.tar.xz

unpack:
tar -xf db-light.tar.xz


mkdir -p bakta_db
tar -xf db-light.tar.xz -C bakta_db


install checkm2 db:
wget https://zenodo.org/record/5571251/files/checkm2_database.tar.gz

mkdir -p checkm2_db
tar -xzf checkm2_database.tar.gz -C checkm2_db

after download there 2 folder. how to move to 1 folder:

mv bakta_db/db-light/* bakta_db/
rmdir bakta_db/db-light


mv checkm2_db/CheckM2_database/* checkm2_db/
rmdir checkm2_db/CheckM2_database