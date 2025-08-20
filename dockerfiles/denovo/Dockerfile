ARG SV_PIPELINE_IMAGE=sv-pipeline:latest
FROM $SV_PIPELINE_IMAGE

ARG DUCKDB_VERSION="1.3.2"
ARG DUCKDB_URI="https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"

COPY /src/denovo /src/denovo/

RUN curl -L -o duckdb_cli-linux-amd64.zip "${DUCKDB_URI}" \
  && unzip duckdb_cli-linux-amd64.zip -d /usr/local/bin \
  && rm duckdb_cli-linux-amd64.zip
