FROM jzj2035198/scfusion:latest

MAINTAINER Zijie Jin (jinzijie@pku.edu.cn)

RUN apt-get update -y \
        && apt-get install -y git

WORKDIR /usr/local/src/

RUN git clone --recursive https://github.com/ZijieJin/scFusion.git

RUN mv scFusion/bin/* scFusion-1.4/bin/ \
        && mv scFusion/scFusion* scFusion-1.4/ \
        && mv scFusion/ scFusion_old/ \
        && mv scFusion-1.4/ scFusion/
