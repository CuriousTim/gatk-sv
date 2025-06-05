ARG SV_PIPELINE_IMAGE=sv-pipeline:latest
FROM $SV_PIPELINE_IMAGE

ARG DUCKDB_VERSION="1.3.0"
ARG DUCKDB_URI="https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"

COPY /src/denovo /src/denovo/

# make & cpp are needed for installing `tidyverse`.
RUN apt-get -qqy update && \
    apt-get -qqy install ca-certificates build-essential && \
    R -e 'install.packages( \
        "https://cran.r-project.org/src/contrib/Archive/yulab.utils/yulab.utils_0.1.9.tar.gz", \
        INSTALL_opts = c("--no-docs", "--without-keep.source"), \
        quiet = TRUE)' && \
    Rscript -e "\
        install.packages( \
            c( \
                'BiocManager', \
                'Rcpp', \
                'tidyverse', \
                'R.utils', \
                'UpSetR', \
                'ggplotify' \
            ), \
            repos = 'https://cran.rstudio.com', \
            quiet = TRUE \
        )" && \
    pip install fsspec gcsfs

RUN curl -L -o duckdb_cli-linux-amd64.zip "${DUCKDB_URI}" \
  && unzip duckdb_cli-linux-amd64.zip -d /usr/local/bin \
  && rm duckdb_cli-linux-amd64.zip
