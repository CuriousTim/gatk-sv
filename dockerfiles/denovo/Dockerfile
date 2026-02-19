ARG DUCKDB_VERSION="1.4.0"
ARG DUCKDB_URI="https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"

ARG BCFTOOLS_VERSION="1.23"
ARG BCFTOOLS_URI="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

ARG HTSLIB_VERSION="1.23"
ARG HTSLIB_URI="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

FROM ubuntu:24.04 AS builder

ARG DUCKDB_VERSION
ARG BCFTOOLS_VERSION
ARG HTSLIB_VERSION
ARG DUCKDB_URI
ARG BCFTOOLS_URI
ARG HTSLIB_URI

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
  --mount=type=cache,target=/var/lib/apt,sharing=locked \
  apt-get update \
  && apt-get install -y --no-install-recommends \
    build-essential \
    libbz2-dev \
    libdeflate-dev \
    liblzma-dev \
    zlib1g-dev \
    unzip

COPY --from=ghcr.io/astral-sh/uv:0.10.4 /uv /uvx /usr/local/bin/
ENV UV_COMPILE_BYTECODE=1
ENV UV_LINK_MODE=copy
ENV UV_NO_DEV=1
ENV UV_PYTHON_INSTALL_DIR=/python
ENV UV_PYTHON_PREFERENCE=only-managed
ENV VIRTUAL_ENV=/opt/venv
RUN uv python install 3.14
ARG PYSAM_VERSION="0.23.3"
ARG PANDAS_VERSION="3.0.1"
RUN uv venv /opt/venv
RUN --mount=type=cache,target=/root/.cache/uv \
  uv pip install "pysam==${PYSAM_VERSION}" "pandas==${PANDAS_VERSION}"

ADD "${DUCKDB_URI}" /tmp/duckdb_cli-linux-amd64.zip
  RUN unzip /tmp/duckdb_cli-linux-amd64.zip -d /usr/local/bin

ADD "${HTSLIB_URI}" /tmp/htslib.tar.bz2
RUN cd /tmp \
  && tar -jxf htslib.tar.bz2 \
  && cd "htslib-${HTSLIB_VERSION}" \
  && ./configure --prefix=/usr/local \
    --with-libdeflate \
  && make \
  && make install

RUN ldconfig

ADD "${BCFTOOLS_URI}" /tmp/bcftools.tar.bz2
RUN cd /tmp \
  && tar -jxf bcftools.tar.bz2 \
  && cd "bcftools-${BCFTOOLS_VERSION}" \
  && ./configure --prefix=/usr/local \
    --with-htslib=/usr/local \
  && make

COPY src/denovo /tmp/denovo

RUN cd /tmp/denovo && make all

FROM ubuntu:24.04

ARG DUCKDB_VERSION
ARG BCFTOOLS_VERSION
ARG HTSLIB_VERSION

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
  --mount=type=cache,target=/var/lib/apt,sharing=locked \
  apt-get update \
  && apt-get install --no-install-recommends -y \
    libbz2-1.0 \
    libdeflate0 \
    liblzma5 \
    zlib1g \
    make

RUN --mount=type=bind,from=builder,src="/tmp/htslib-${HTSLIB_VERSION}",target=/tmp/htslib \
  --mount=type=bind,from=builder,src="/tmp/bcftools-${BCFTOOLS_VERSION}",target=/tmp/bcftools \
  --mount=type=bind,from=builder,src=/tmp/denovo,target=/tmp/denovo \
  cd /tmp/htslib && make install \
  && cd /tmp/bcftools && make install \
  && cd /tmp/denovo && make install

COPY --from=builder /usr/local/bin/duckdb /usr/local/bin/duckdb
COPY --from=builder --chown=python:python /python /python
COPY --from=builder /opt/venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

COPY src/denovo/null_low_coverage_gt.py /opt/gatk-sv/denovo/null_low_coverage_gt.py

RUN ldconfig

RUN groupadd --system --gid 999 nonroot \
 && useradd --system --gid 999 --uid 999 --create-home nonroot
USER nonroot
WORKDIR /tmp

CMD ["bash"]
