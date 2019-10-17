FROM nfcore/base

MAINTAINER Diego Brambilla
LABEL description="Docker image containing all requirements for nf-core/ampliseq pipeline, including dependecies for R-DADA2 scripts"

## Install necessary tools
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a && \
    cd /opt && git clone https://github.com/erikrikarddaniel/eemisdada2.git && \
    ln -s /opt/eemisdada2/src/R/dada2* /usr/local/bin
RUN mkdir -p /data/taxonomy && \
    wget https://files.plutof.ut.ee/public/orig/D6/96/D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip && \
    unzip D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip sh_general_release_dynamic_all_02.02.2019.fasta.gz && \
    gzip sh_general_release_dynamic_all_02.02.2019.fasta.gz && \
    rm D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip && \
    ln -s /opt/sh_general_release_dynamic_all_02.02.2019.fasta.gz /data/taxonomy/UNITE.fna
RUN wget -O UNITE_v2019_July2019.RData http://www2.decipher.codes/Classification/TrainingSets/UNITE_v2019_July2019.RData && \
    ln -s  /opt/UNITE_v2019_July2019.RData /data/taxonomy/UNITE.RData
RUN wget -O silva_nr_v132_train_set.fa.gz https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1 &&\
   ln -s /opt/silva_nr_v132_train_set.fa.gz /data/taxonomy/SILVA.fna
RUN wget -O silva_species_assignment_v132.fa.gz https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1 && \
    ln -s /opt/silva_species_assignment_v132.fa.gz /data/taxonomy/SILVA_species.fna
RUN wget -O SILVA_SSU_r132_March2018.RData http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r132_March2018.RData && \
    ln -s /opt/SILVA_SSU_r132_March2018.RData /data/taxonomy/SILVA.RData
RUN wget -O GTDB_bac-arc_ssu_r86.fa.gz https://zenodo.org/record/2658728/files/GTDB_bac-arc_ssu_r86.fa.gz?download=1 &&\
   ln -s /opt/GTDB_bac-arc_ssu_r86.fa.gz /data/taxonomy/GTDB.fna
RUN wget -O GTDB_dada2_assignment_species.fa.gz https://zenodo.org/record/2658728/files/GTDB_dada2_assignment_species.fa.gz?download=1 && \
   ln -s /opt/GTDB_dada2_assignment_species.fa.gz /data/taxonomy/GTDB_species.fna
RUN wget -O GTDB_r89-mod_June2019.RData http://www2.decipher.codes/Classification/TrainingSets/GTDB_r89-mod_June2019.RData && \
   ln -s /opt/GTDB_r89-mod_June2019.RData /data/taxonomy/GTDB.RData
RUN wget -O pr2_version_4.12.0_18S_dada2.fasta.gz https://github.com/pr2database/pr2database/releases/download/v4.12.0/pr2_version_4.12.0_18S_dada2.fasta.gz
    ln -s /opt/pr2_version_4.12.0_18S_dada2.fasta.gz /data/taxonomy/PR2.fna
ENV PATH /opt/conda/envs/bioatlas-ampliflow/bin:$PATH

