FROM nfcore/base

MAINTAINER Diego Brambilla
LABEL description="Docker image containing all requirements for nf-core/ampliseq pipeline, including dependecies for R-DADA2 scripts"

## Install necessary tools
COPY environment.yml /
RUN apt-get update && \
    apt-get install -y --no-install-recommends unzip \
    && apt-get clean
RUN conda env create -f /environment.yml && conda clean -a && \
    cd /opt && git clone https://github.com/erikrikarddaniel/eemisdada2.git && \
    ln -s /opt/eemisdada2/src/R/dada2* /usr/local/bin
RUN mkdir -p /data/taxonomy && cd /opt &&\
    wget https://files.plutof.ut.ee/public/orig/EB/0C/EB0CCB3A871B77EA75E472D13926271076904A588D2E1C1EA5AFCF7397D48378.zip && \
    unzip -p "EB0CCB3A871B77EA75E472D13926271076904A588D2E1C1EA5AFCF7397D48378.zip" > "sh_general_release_02.02.2019.fasta.gz" && \
    gzip sh_general_release_02.02.2019.fasta.gz && \
    rm EB0CCB3A871B77EA75E472D13926271076904A588D2E1C1EA5AFCF7397D48378.zip && \
    ln -s /opt/sh_general_release_02.02.2019.fasta.gz /data/taxonomy/UNITE.fna.gz
RUN cd /opt && wget -O UNITE_v2019_July2019.RData http://www2.decipher.codes/Classification/TrainingSets/UNITE_v2019_July2019.RData && \
    ln -s  /opt/UNITE_v2019_July2019.RData /data/taxonomy/UNITE.RData
RUN cd /opt && wget -O silva_nr_v132_train_set.fa.gz https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1 &&\
   ln -s /opt/silva_nr_v132_train_set.fa.gz /data/taxonomy/SILVA.fna.gz
RUN cd /opt && wget -O silva_species_assignment_v132.fa.gz https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1 && \
    ln -s /opt/silva_species_assignment_v132.fa.gz /data/taxonomy/SILVA_species.fna.gz
RUN cd /opt && wget -O SILVA_SSU_r132_March2018.RData http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r132_March2018.RData && \
    ln -s /opt/SILVA_SSU_r132_March2018.RData /data/taxonomy/SILVA.RData
RUN cd /opt && wget -O GTDB_bac-arc_ssu_r86.fa.gz https://zenodo.org/record/2541239/files/GTDB_bac-arc_ssu_r86.fa.gz?download=1 &&\
   ln -s /opt/GTDB_bac-arc_ssu_r86.fa.gz /data/taxonomy/GTDB.fna.gz
RUN cd /opt && wget -O GTDB_dada2_assignment_species.fa.gz https://zenodo.org/record/2658728/files/GTDB_dada2_assignment_species.fa.gz?download=1 && \
   ln -s /opt/GTDB_dada2_assignment_species.fa.gz /data/taxonomy/GTDB_species.fna.gz
RUN cd /opt && wget -O GTDB_r89-mod_June2019.RData http://www2.decipher.codes/Classification/TrainingSets/GTDB_r89-mod_June2019.RData && \
   ln -s /opt/GTDB_r89-mod_June2019.RData /data/taxonomy/GTDB.RData
RUN cd /opt && wget -O pr2_version_4.12.0_18S_dada2.fasta.gz https://github.com/pr2database/pr2database/releases/download/v4.12.0/pr2_version_4.12.0_18S_dada2.fasta.gz && \
   ln -s /opt/pr2_version_4.12.0_18S_dada2.fasta.gz /data/taxonomy/PR2.fna.gz
ENV PATH /opt/conda/envs/bioatlas-ampliflow/bin:$PATH

