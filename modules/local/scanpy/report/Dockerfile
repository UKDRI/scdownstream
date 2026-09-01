# Image: scanpy_1.11.4_coreinf_0.3
#
FROM gcfntnu/scanpy:1.11.4

RUN apt update

RUN wget -qO /tmp/quarto-1.9.38-linux-amd64.deb https://github.com/quarto-dev/quarto-cli/releases/download/v1.9.38/quarto-1.9.38-linux-amd64.deb && \
    apt install -y /tmp/quarto-1.9.38-linux-amd64.deb && \
    rm /tmp/quarto-1.9.38-linux-amd64.deb && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir papermill
RUN pip install --no-cache-dir nbconvert
RUN pip install --no-cache-dir gprofiler-official
RUN pip install --no-cache-dir liana==1.6.1