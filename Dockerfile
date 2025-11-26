FROM mambaorg/micromamba:2.4.0

# Create environment
# conda env export --from-history -n arvia > env.yaml # the "--from-history" is important, only packages i asked for
COPY --chown=$MAMBA_USER:$MAMBA_USER arvia/envs/arvia.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Copy code into docker
# COPY . .

# Flag to activate conda inside the build command
# (this is set automatically with "docker run" but not "docker build"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Test tool and update db
RUN arvia --help && \
    arvia dbs -o ./arvia_dbs && \
    arvia test -o ./arvia_test && \
    rm -r ./arvia_dbs ./arvia_test

# Sanity check
RUN arvia --help

 # /opt/conda/share/amrfinderplus/data/latest # amrfinder
 # /opt/conda/db # mlst

