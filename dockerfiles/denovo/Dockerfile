FROM ubuntu:24.04

ARG DUCKDB_VERSION="1.4.0"
ARG DUCKDB_URI="https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"

ARG BCFTOOLS_VERSION="1.22"
ARG BCFTOOLS_URI="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

ARG HTSLIB_VERSION="1.22.1"
ARG HTSLIB_URI="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    build-essential \
    libbz2-dev \
    libdeflate-dev \
    liblzma-dev \
    zlib1g-dev \
    unzip \
  && rm -rf /var/lib/apt/lists/*

ADD "${DUCKDB_URI}" /tmp/duckdb_cli-linux-amd64.zip
RUN unzip /tmp/duckdb_cli-linux-amd64.zip -d /usr/local/bin \
  && rm /tmp/duckdb_cli-linux-amd64.zip

ADD "${HTSLIB_URI}" /tmp/htslib.tar.bz2
RUN cd /tmp \
  && tar -jxf htslib.tar.bz2 \
  && cd "htslib-${HTSLIB_VERSION}" \
  && ./configure --prefix=/usr/local \
    --with-libdeflate \
  && make \
  && make install \
  && cd .. \
  && rm -rf "htslib-${HTSLIB_VERSION}" htslib.tar.bz2

RUN ldconfig

ADD "${BCFTOOLS_URI}" /tmp/bcftools.tar.bz2
RUN cd /tmp \
  && tar -jxf bcftools.tar.bz2 \
  && cd "bcftools-${BCFTOOLS_VERSION}" \
  && ./configure --prefix=/usr/local \
    --with-htslib=/usr/local \
  && make \
  && make install \
  && cd .. \
  && rm -rf "bcftools-${BCFTOOLS_VERSION}" bcftools.tar.bz2

RUN --mount=type=bind,source=src/denovo,target=/tmp/denovo,rw \
  cd /tmp/denovo \
  && make all \
  && make install

CMD ["bash"]
