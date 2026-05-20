FROM mambaorg/micromamba:1.5.0

LABEL maintainer="SNPhylo2 Development Team"
LABEL description="SNPhylo2: Next-Generation Phylogenomic Pipeline"
LABEL version="0.1.0"

# Switch to root for installation
USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install external bioinformatics tools
# IQ-TREE2
RUN wget -q https://github.com/iqtree/iqtree2/releases/download/v2.3.4/iqtree-2.3.4-Linux.tar.gz \
    && tar -xzf iqtree-2.3.4-Linux.tar.gz \
    && mv iqtree-2.3.4-Linux/bin/iqtree2 /usr/local/bin/ \
    && rm -rf iqtree-2.3.4-Linux*

# RAxML-NG
RUN wget -q https://github.com/amkozlov/raxml-ng/releases/download/1.2.0/raxml-ng_v1.2.0_linux_x86_64.zip \
    && unzip -q raxml-ng_v1.2.0_linux_x86_64.zip \
    && mv raxml-ng /usr/local/bin/ \
    && rm raxml-ng_v1.2.0_linux_x86_64.zip

# FastTree
RUN wget -q http://www.microbesonline.org/fasttree/FastTree \
    && chmod +x FastTree \
    && mv FastTree /usr/local/bin/

# PLINK2
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240123.zip \
    && unzip -q plink2_linux_avx2_20240123.zip \
    && mv plink2 /usr/local/bin/ \
    && rm plink2_linux_avx2_20240123.zip

# BCFtools
RUN wget -q https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
    && tar -xjf bcftools-1.19.tar.bz2 \
    && cd bcftools-1.19 \
    && ./configure --prefix=/usr/local \
    && make && make install \
    && cd .. && rm -rf bcftools-1.19*

# HTSlib (for tabix)
RUN wget -q https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 \
    && tar -xjf htslib-1.19.tar.bz2 \
    && cd htslib-1.19 \
    && ./configure --prefix=/usr/local \
    && make && make install \
    && cd .. && rm -rf htslib-1.19*

# Set working directory
WORKDIR /app

# Copy conda environment file
COPY environment.yml /app/

# Create conda environment
RUN micromamba install -y -n base -f /app/environment.yml \
    && micromamba clean --all --yes

# Copy SNPhylo2 source code
COPY snphylo2/ /app/snphylo2/
COPY pyproject.toml /app/
COPY README.md /app/
COPY LICENSE /app/

# Install SNPhylo2
RUN pip install -e . --no-deps

# Create directories for data and results
RUN mkdir -p /data /results

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV MPLBACKEND=Agg

# Switch back to micromamba user
USER $MAMBA_USER

# Default command
ENTRYPOINT ["snphylo2"]
CMD ["--help"]
