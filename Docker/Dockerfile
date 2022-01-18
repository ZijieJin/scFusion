FROM jzj2035198/scfusion:latest

MAINTAINER Zijie Jin (jinzijie@pku.edu.cn)

RUN apt-get update -y \
        && apt-get install -y git

WORKDIR /usr/local/src/

RUN git clone --recursive https://github.com/ZijieJin/scFusion.git

RUN unzip scFusion/data/hg19mappability75.txt.zip -d scFusion/data/

RUN mv scFusion-1.4/data/hg19.fa scFusion/data \
        && mv scFusion-1.4/data/ref_annot.gtf scFusion/data/
