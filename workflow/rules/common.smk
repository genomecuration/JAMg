import os

from snakemake.utils import validate

validate(config, "../config/config.schema.yaml")

GENOME = config["genome"]
GENOME_BASENAME = os.path.basename(GENOME)
OUTDIR = config["outdir"]
THREADS = config.get("threads", 8)


def pasa_tmp_dir():
    return f"/dev/shm/{os.environ.get('USER', 'jamg')}"
