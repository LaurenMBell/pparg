FROM mambaorg/micromamba:1.5.8

WORKDIR /app

# Copy environment spec first for better layer caching
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /app/environment.yml

# Create conda env (includes RDKit)
RUN micromamba create -y -n pparg -f /app/environment.yml && \
    micromamba clean -a -y

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV MPLCONFIGDIR=/tmp/mplconfig
ENV PYTHONUNBUFFERED=1

# Copy app code
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

EXPOSE 10000

# Render sets $PORT. Gunicorn binds to it.
CMD micromamba run -n pparg gunicorn -w 2 -k gthread --threads 4 -b 0.0.0.0:${PORT:-10000} app:app

