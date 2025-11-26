# gast_project


### Prerequisites

* Install required packages for automated data retrieval.
* Download the associated databases (Checkm2 and bakta).


---

**Needed reference databases installation**
These needs to be installed before workflow running

install bakta_db:
wget https://zenodo.org/record/14916843/files/db-light.tar.xz

# unpack
tar -xf db-light.tar.xz

# rename to bakta_db
mv db-light bakta_db


install checkm2 db:
wget https://zenodo.org/record/5571251/files/checkm2_database.tar.gz

mkdir -p /path/to/checkm2_db
tar -xzf checkm2_database.tar.gz -C /path/to/checkm2_db