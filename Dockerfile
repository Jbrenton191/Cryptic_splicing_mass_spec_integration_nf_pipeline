# -------------------------------------------------------------------
# dasper (Bioconductor 3.17) in an isolated conda env, alongside
# rocker's R 4.6.1, which stays the default R / Rscript.
# -------------------------------------------------------------------
FROM rocker/r-ver:4.6.1

# -------------------------------------------------------------------
# 0. Build settings
# -------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# -------------------------------------------------------------------
# 1. System dependencies
# -------------------------------------------------------------------
RUN apt-get update && apt-get install -y \
      libpng-dev \
      libjpeg-dev \
      libtiff-dev \
      libcairo2-dev \
      libglib2.0-dev \
      libfreetype6-dev \
      libx11-dev \
      libxt-dev \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      build-essential \
      libgit2-dev \
      libgsl-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      pkg-config \
      git \
      vim \
      wget \
      tmux \
      rename \
    && rm -rf /var/lib/apt/lists/*

# -------------------------------------------------------------------
# 2. Install Miniforge (conda + mamba)
#    No `conda init`, base env never auto-activates, and conda goes at
#    the END of PATH, so rocker's R is always found first.
# -------------------------------------------------------------------
ENV CONDA_DIR=/opt/conda
RUN wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" \
      -O /tmp/miniforge.sh \
 && bash /tmp/miniforge.sh -b -p "${CONDA_DIR}" \
 && rm /tmp/miniforge.sh \
 && "${CONDA_DIR}/bin/conda" config --system --set always_yes yes --set changeps1 yes --set auto_activate_base false \
 && "${CONDA_DIR}/bin/conda" clean -afy
ENV PATH="${PATH}:${CONDA_DIR}/bin"

# -------------------------------------------------------------------
# 3. Create the conda env with R 4.3.0
# -------------------------------------------------------------------
RUN mamba create -n dasper_env -c conda-forge -c bioconda \
      conda-forge::r-base=4.3.0 \
 && mamba clean -afy

# -------------------------------------------------------------------
# 4. Install BiocManager and set Bioconductor 3.17 (in the conda R)
#    Called by full path; R_HOME (set by rocker) is unset for this
#    command only so conda's R uses its own installation.
# -------------------------------------------------------------------
RUN env -u R_HOME /opt/conda/envs/dasper_env/bin/Rscript \
      -e "if (!require(BiocManager)) install.packages('BiocManager', repos='https://www.stats.bris.ac.uk/R/')" \
      -e "BiocManager::install(version = '3.17', force = TRUE, ask = FALSE)"

# -------------------------------------------------------------------
# 5. Install dependencies from conda-forge / bioconda
# -------------------------------------------------------------------
RUN mamba install -n dasper_env -c conda-forge -c bioconda \
      conda-forge::r-recommended \
      conda-forge::r-ggpubr \
      bioconda::bioconductor-basilisk \
      bioconda::bioconductor-biocfilecache \
      bioconda::bioconductor-genomeinfodb \
      bioconda::bioconductor-genomicfeatures \
      bioconda::bioconductor-genomicranges \
      bioconda::bioconductor-megadepth \
      bioconda::bioconductor-plyranges \
      conda-forge::r-reticulate \
      bioconda::bioconductor-rtracklayer \
      conda-forge::r-ragg \
      --freeze-installed \
 && mamba clean -afy

# -------------------------------------------------------------------
# 6. Install dasper and tidyverse (in the conda R)
# -------------------------------------------------------------------
RUN env -u R_HOME /opt/conda/envs/dasper_env/bin/Rscript \
      -e "if (!requireNamespace('dasper', quietly = TRUE)) BiocManager::install('dasper', ask = FALSE, update = FALSE)" \
      -e "if (!requireNamespace('tidyverse', quietly = TRUE)) install.packages('tidyverse', repos = 'https://www.stats.bris.ac.uk/R/')" \
      -e "if (!requireNamespace('optparse', quietly = TRUE)) install.packages('optparse', repos = 'https://www.stats.bris.ac.uk/R/')"

# -------------------------------------------------------------------
# 7. Check both R installations (build fails if either is wrong)
# -------------------------------------------------------------------
RUN test "$(command -v R)" = "/usr/local/bin/R" \
 && test "$(command -v Rscript)" = "/usr/local/bin/Rscript" \
 && Rscript -e 'cat("default R:", R.version.string, "\n")' \
 && env -u R_HOME /opt/conda/envs/dasper_env/bin/Rscript \
      -e 'library(dasper); cat("dasper", as.character(packageVersion("dasper")), "on", R.version.string, "\n")'
# -------------------------------------------------------------------
# 8. Pipeline code (cloned into the home folder)
# -------------------------------------------------------------------
RUN git clone https://github.com/Jbrenton191/Cryptic_splicing_mass_spec_integration_nf_pipeline.git \
      /root/Cryptic_splicing_mass_spec_integration_nf_pipeline

# -------------------------------------------------------------------
# 9. Nextflow (own env; brings its own Java)
# -------------------------------------------------------------------
RUN mamba create -n nextflow_env -c conda-forge -c bioconda nextflow \
 && mamba clean -afy
ENV PATH="${PATH}:${CONDA_DIR}/envs/nextflow_env/bin"

# -------------------------------------------------------------------
# 10. Start settings: activate nextflow_env, start in the pipeline folder, run bash
# -------------------------------------------------------------------
RUN echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate nextflow_env" >> /root/.bashrc



