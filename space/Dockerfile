# AURORA gene expression viewer.
#
# The viewer decodes gene-expression fields from the original prediction HDF5
# files on demand, so no model weights or pre-rendered images are baked in.
# Point it at a Hugging Face dataset with AURORA_BACKEND=hf, or at a folder of
# files with AURORA_BACKEND=local / AURORA_LOCAL_DIR.

FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# imagecodecs ships the zstd codec used to decompress the prediction chunks, and
# h5py/numpy/pillow all publish manylinux wheels, so no compiler is needed.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Hugging Face Spaces routes to app_port; anything that binds 0.0.0.0 works.
# The LUAD predictions are split across two folders in the dataset repo (112 um
# in TCGA_LUAD_website, 7/28 um in predictions_7um_28um/TCGA_LUAD), so the
# subfolder list is comma-separated.  The H&E pyramids live in HE_image_website.
ENV AURORA_BACKEND=hf \
    AURORA_HF_REPO=AURORAData/AURORA_prediction \
    AURORA_HF_SUBFOLDER=TCGA_LUAD_website,predictions_7um_28um/TCGA_LUAD \
    AURORA_HE_REPO=AURORAData/AURORA_prediction \
    AURORA_HE_SUBFOLDER=HE_image_website \
    AURORA_DEFAULT_SAMPLE=TCGA-44-5645-01B-04-BS4 \
    AURORA_DEFAULT_GENE=EPCAM \
    AURORA_PREFETCH=1 \
    PORT=7860

EXPOSE 7860

# One worker: the decode cache is per-process, and a second worker would double
# the memory for no benefit at this traffic level.
CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-7860} --workers 1 --timeout-keep-alive 30"]
