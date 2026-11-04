# Solver of the IJES revision: Marmot and EdelweissFE at the tags the analyses were run with, the plug-ins of
# Marmot/plugins/ijes-revision, and the conda environment of the workstation (env/linux-64.explicit.txt).
#   docker build -t ijes-interface .        or    apptainer pull docker://ghcr.io/alexsta-collab/marmot:ijes-reproducibility-2026-11
FROM mambaorg/micromamba:2.0.5

ARG MARMOT_REPO=https://github.com/AlexSTA-collab/Marmot
ARG EDELWEISSFE_REPO=https://github.com/AlexSTA-collab/EdelweissFE
ARG MARMOT_TAG=ijes-revision-2026-11c
ARG EDELWEISSFE_TAG=ijes-revision-2026-11

USER root
ENV ENV_PREFIX=/opt/env MARMOT_PREFIX=/opt/marmot PATH=/opt/env/bin:$PATH

# the environment, package by package as on the workstation
COPY env/linux-64.explicit.txt /tmp/linux-64.explicit.txt
RUN micromamba create -y -p $ENV_PREFIX --file /tmp/linux-64.explicit.txt && micromamba clean -a -y

# Marmot, as ~/psie26/build_marmot.sh built it
ADD ${MARMOT_REPO}/archive/refs/tags/${MARMOT_TAG}.tar.gz /tmp/marmot.tar.gz
RUN mkdir -p /opt/src/Marmot && tar -xzf /tmp/marmot.tar.gz -C /opt/src/Marmot --strip-components=1 \
 && cmake -S /opt/src/Marmot -B /tmp/marmot-build -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=$ENV_PREFIX/bin/c++ -DCMAKE_C_COMPILER=$ENV_PREFIX/bin/cc -DCMAKE_PREFIX_PATH=$ENV_PREFIX \
      -DCMAKE_INSTALL_PREFIX=$MARMOT_PREFIX -DBUILD_TESTING=OFF -DMARMOT_BUILD_PYTHON_BINDINGS=OFF \
 && cmake --build /tmp/marmot-build -j "$(nproc)" && cmake --install /tmp/marmot-build \
 && CXX=$ENV_PREFIX/bin/c++ /opt/src/Marmot/plugins/ijes-revision/build.sh $MARMOT_PREFIX $ENV_PREFIX $MARMOT_PREFIX/plugins \
 && rm -rf /tmp/marmot-build /tmp/marmot.tar.gz

# EdelweissFE, linked to that Marmot
ADD ${EDELWEISSFE_REPO}/archive/refs/tags/${EDELWEISSFE_TAG}.tar.gz /tmp/edelweissfe.tar.gz
RUN mkdir -p /opt/src/EdelweissFE && tar -xzf /tmp/edelweissfe.tar.gz -C /opt/src/EdelweissFE --strip-components=1 \
 && cd /opt/src/EdelweissFE && MARMOT_INSTALL_DIR=$MARMOT_PREFIX $ENV_PREFIX/bin/python -m pip install --no-deps --no-build-isolation . \
 && rm -rf /tmp/edelweissfe.tar.gz

# the package: job files, scorer, figure scripts
COPY . /opt/ijes
RUN $ENV_PREFIX/bin/python -m pip install --no-deps --no-build-isolation -e /opt/ijes && chmod +x /opt/ijes/docker/ijes-run
ENV PATH=/opt/ijes/docker:$PATH OMP_NUM_THREADS=4 PYTHONUNBUFFERED=1
WORKDIR /work
CMD ["ijes", "--help"]
