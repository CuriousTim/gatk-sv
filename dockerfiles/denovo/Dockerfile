FROM ubuntu:24.04

ARG DUCKDB_VERSION="1.3.2"
ARG DUCKDB_URI="https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"

# COPY /src/denovo /src/denovo/

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    unzip \
  && rm -rf /var/lib/apt/lists/*

ADD $DUCKDB_URI /tmp/duckdb_cli-linux-amd64.zip
RUN unzip /tmp/duckdb_cli-linux-amd64.zip -d /usr/local/bin \
  && rm /tmp/duckdb_cli-linux-amd64.zip
