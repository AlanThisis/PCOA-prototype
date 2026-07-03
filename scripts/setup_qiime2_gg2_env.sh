#!/usr/bin/env bash
set -euo pipefail

QIIME2_VERSION="${QIIME2_VERSION:-2026.4}"
ENV_NAME="${QIIME2_ENV_NAME:-rachis-qiime2-${QIIME2_VERSION}}"
BASE_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/${QIIME2_VERSION}/qiime2/released"

case "$(uname -s)" in
  Linux)
    ENV_FILE_URL="${BASE_URL}/rachis-qiime2-linux-64-conda.yml"
    conda env create --name "${ENV_NAME}" --file "${ENV_FILE_URL}"
    ;;
  Darwin)
    ENV_FILE_URL="${BASE_URL}/rachis-qiime2-osx-64-conda.yml"
    CONDA_SUBDIR=osx-64 conda env create --name "${ENV_NAME}" --file "${ENV_FILE_URL}"

    CONDA_BASE="$(conda info --base)"
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}"
    conda config --env --set subdir osx-64
    conda deactivate
    ;;
  *)
    echo "Unsupported OS: $(uname -s). Use Linux/WSL or macOS." >&2
    exit 1
    ;;
esac

conda env update --name "${ENV_NAME}" --file environment.yml
conda run --name "${ENV_NAME}" python -m pip install --force-reinstall --no-deps setuptools==80.9.0
conda run --name "${ENV_NAME}" python -c "import pkg_resources; print('pkg_resources OK')"
conda run --name "${ENV_NAME}" python -m pip install fastq-dl==4.0.1
conda run --name "${ENV_NAME}" python -m pip install --no-deps redbiom==0.3.9 q2-greengenes2==2024.1

conda run --name "${ENV_NAME}" qiime info
conda run --name "${ENV_NAME}" qiime greengenes2 --help >/dev/null
conda run --name "${ENV_NAME}" fastq-dl --version
conda run --name "${ENV_NAME}" seqkit version
conda run --name "${ENV_NAME}" deblur --version
conda run --name "${ENV_NAME}" python -c "import biom, skbio, pandas, matplotlib, requests; print('Python deps OK')"

cat <<EOF

Environment ready:
  conda activate ${ENV_NAME}

GG2 reference artifacts are not bundled. Download them with the README command
before running src/unifrac.py.
EOF
