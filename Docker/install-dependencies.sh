#!/usr/bin/env bash
export BWA_VERSION=${BWA_VERSION:-0.7.18}
export JUICER_TOOLS_VERSION=${JUICER_TOOLS_VERSION:-2.20.00}
export PAIRIX_VERSION=${PAIRIX_VERSION:-0.3.6}
export SAMTOOLS_VERSION=${SAMTOOLS_VERSION:-1.21}

# For sorting, LC_ALL is C
export LC_ALL=C

apt update && apt install -y \
    bzip2 curl gawk gcc make git wget aria2 mc unzip sudo \
    libbz2-dev libz-dev locales libcurl4-openssl-dev \
    openjdk-8-jdk python3 python-is-python3 python3-pip \
&& apt clean && rm -rf /var/lib/apt/list/* \
&& echo 'kasm-user ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers 

echo 'alias awk=gawk' >> /etc/profile.d/aliases.sh && chmod +x /etc/profile.d/aliases.sh && \
    locale-gen en_US.UTF-8 

cd /opt/ || exit

# Install BWA
curl -OL "https://github.com/lh3/bwa/archive/v${BWA_VERSION}.zip" && \
    unzip "v${BWA_VERSION}.zip" && \
    cd "bwa-${BWA_VERSION}/" && \
    make && \
    cp bwa /usr/local/bin && \
    cd .. && \
    rm -rf "bwa-v${BWA_VERSION}"

# Install Samtools
curl -OL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    bunzip2 "samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar xf "samtools-${SAMTOOLS_VERSION}.tar" && \
    cd "samtools-${SAMTOOLS_VERSION}" && \
    ./configure --without-curses --disable-bz2 --disable-lzma && \
    make && \
    make install && \
    cd .. && \
    rm -rf "samtools-${SAMTOOLS_VERSION}" "samtools-${SAMTOOLS_VERSION}.tar"

# Install pairix, only need bam2pairs + deps
curl -OL "https://github.com/4dn-dcic/pairix/archive/${PAIRIX_VERSION}.zip" && \
    unzip "${PAIRIX_VERSION}.zip" && \
    cd "pairix-${PAIRIX_VERSION}/" && \
    make && \
    cp bin/pairix bin/bgzip /usr/local/bin && \
    cp util/bam2pairs/bam2pairs /opt && \
    cd .. && \
    rm -rf "${PAIRIX_VERSION}.zip" "pairix-${PAIRIX_VERSION}"

# Install Juicer
git clone --branch encode https://github.com/theaidenlab/juicer.git && \    
    chmod +x juicer/CPU/* juicer/CPU/common/* juicer/misc/* 
#    && \
#    find -mindepth 1 -maxdepth 1  -type d -not -name "CPU" -not -name ".git" -not -name "misc" | xargs rm -rf

# Install Juicer tools
wget "https://github.com/aidenlab/Juicebox/releases/download/v${JUICER_TOOLS_VERSION}/juicer_tools.${JUICER_TOOLS_VERSION}.jar" -O /opt/juicer/CPU/common/juicer_tools.jar && \
    chmod 666 /opt/juicer/CPU/common/juicer_tools.jar && \
    ln -s juicer/CPU scripts && \
    wget "https://s3.us-east-1.wasabisys.com/hicfiles/public/Juicebox/Juicebox-2.15.jar" -O /opt/juicer/CPU/common/Juicebox.jar && \
    chmod a+x /opt/juicer/CPU/common/Juicebox.jar && \
    if [ -d "$HOME/Desktop" ]; then \
        ln -s /opt/juicer/CPU/common/Juicebox.jar $HOME/Desktop/Juicebox.jar && \
        ln -s /opt/juicer/CPU/common/juicer_tools.jar $HOME/Desktop/juicer_tools.jar; \
    fi && \
    mkdir -p /aidenlab && chown -R 1000:0 /aidenlab && chmod -R a+rwX /aidenlab && ln -s /opt/juicer/CPU /aidenlab/scripts 

# Install latest VS Code cli
curl -fsSL "https://code.visualstudio.com/sha/download?build=stable&os=cli-alpine-x64" -o vscode_cli_alpine_x64.tar.gz && \
    tar xzf vscode_cli_alpine_x64.tar.gz && \
    mv code /usr/local/bin && \
    rm vscode_cli_alpine_x64.tar.gz 

# Install hic-straw, pandas, numpy and toolshed python packages
pip install pandas numpy toolshed pycurl hic-straw 
