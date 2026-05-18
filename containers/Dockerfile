# JAMg v2 container image -- Docker equivalent of containers/jamg.def +
# containers/jamg-base.def collapsed into a single Dockerfile (Docker's own
# layer cache handles the equivalent of the two-SIF split).
#
# Build from repo root:
#   docker build -f containers/Dockerfile -t jamg:2.0.0a1 .

# Base layer (mirrors containers/jamg-base.def): debian-slim + apt bootstrap
# tools + statically linked micromamba.
FROM debian:trixie-slim
# Was debian:bookworm-slim; bumped to trixie (glibc 2.40) so the SIF satisfies
# apptainer's libfakeroot.so GLIBC_2.38 requirement -- bookworm's glibc 2.36
# was too old. Same constraint applies to Docker builds running under apptainer
# inspect or rootless docker with apparmor.

# %post equivalent (base layer)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
         ca-certificates curl bzip2 tar python3 \
    && rm -rf /var/lib/apt/lists/* \
    && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
         | tar -xj -C /usr/local bin/micromamba \
    && chmod +x /usr/local/bin/micromamba \
    && mkdir -p /opt/conda

ENV MAMBA_ROOT_PREFIX=/opt/conda

# %files equivalent (main layer)
COPY containers/environment.yml /opt/environment.yml
COPY containers/post-install.sh /opt/post-install.sh
COPY bin/                        /opt/jamg/bin/
COPY PerlLib/                    /opt/jamg/PerlLib/

# %post equivalent (main layer) -- mirrors jamg.def line for line.
# Reviewer C1/C2/C3 fixes: export MAMBA_ROOT_PREFIX explicitly, PATH set, activate base.
# Reviewer I4: clean AFTER post-install.
RUN micromamba install -n base -y -f /opt/environment.yml \
    && export PATH=/opt/conda/bin:$PATH \
    && eval "$(micromamba shell hook -s bash)" \
    && micromamba activate base \
    && chmod +x /opt/post-install.sh \
    && bash /opt/post-install.sh \
    && micromamba clean -ay \
    && echo "export PERL5LIB=/opt/jamg/PerlLib" >> /etc/profile.d/jamg.sh \
    && echo "export PATH=/opt/jamg/bin:\$PATH"   >> /etc/profile.d/jamg.sh

# %environment equivalent -- PATH order /opt/jamg/bin first (reviewer I3).
ENV PERL5LIB=/opt/jamg/PerlLib
ENV PATH=/opt/jamg/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV LC_ALL=C.UTF-8

# ENTRYPOINT intentionally OMITTED in Phase 1 (matches jamg.def -- bin/jamg
# is a Phase 2 deliverable). Code-review finding B1: pointing the entrypoint
# at a non-existent file would make `docker run jamg:2.0.0a1` fail on every
# invocation. Phase 2 will re-add:
#     ENTRYPOINT ["/opt/jamg/bin/jamg"]
#     CMD []
