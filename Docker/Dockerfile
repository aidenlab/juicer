#FROM ubuntu:jammy
FROM nvidia/cuda:11.7.1-devel-ubuntu22.04

ENV BWA_VERSION=0.7.18
ENV JUICER_TOOLS_VERSION=2.20.00
ENV PAIRIX_VERSION=0.3.6
ENV SAMTOOLS_VERSION=1.21

# For sorting, LC_ALL is C
ENV LC_ALL=C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:/opt/juicer/misc:/aidenlab:$PATH
WORKDIR /opt/
COPY install-dependencies.sh /opt/install-dependencies.sh
# Copy scripts for downloading demo data
COPY download-and-run-demo.sh /aidenlab/
COPY download-demo.txt /aidenlab/
RUN chmod +x /opt/install-dependencies.sh && /opt/install-dependencies.sh && \
    chmod +x /aidenlab/download-and-run-demo.sh && \
    ln -s /aidenlab/download-and-run-demo.sh /usr/bin/download-and-run-demo 
ENTRYPOINT [ "/opt/scripts/juicer.sh" ]