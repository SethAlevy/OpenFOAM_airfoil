FROM ubuntu:22.04

USER root
ENV DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC
ENV POETRY_VIRTUALENVS_CREATE=true
ENV PATH="/root/.local/bin:$PATH"

RUN apt-get update && \
    apt-get install -y \
      python3.11 python3.11-venv python3-pip \
      git curl wget unzip build-essential ca-certificates \
      jq dos2unix && \
    rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 && \
    python3 -m pip install --upgrade pip && \
    python3 -m pip install poetry

# Install OpenFOAM v2406 and build cfMesh BEFORE copying app (keeps cache)
RUN wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | bash && \
    apt-get update && apt-get install -y openfoam2406 openfoam2406-dev && \
    ls -la /usr/lib/openfoam/openfoam2406/wmake/scripts || (echo "wmake scripts missing" && exit 1)

RUN curl -L --fail -H "User-Agent: curl" \
      https://develop.openfoam.com/Community/integration-cfmesh/-/archive/v2406/integration-cfmesh-v2406.zip \
      -o /tmp/cfmesh-install.zip

RUN unzip -q /tmp/cfmesh-install.zip -d /opt && \
    mv /opt/integration-cfmesh-v2406 /opt/cfmesh-src && \
    find /opt/cfmesh-src -type f -name "Allwmake*" -exec chmod +x {} \;

RUN /usr/bin/openfoam2406 -c 'cd /opt/cfmesh-src && /bin/bash ./Allwmake'

WORKDIR /app

# Cache-friendly deps
COPY pyproject.toml poetry.lock* /app/
RUN mkdir -p /root/.cache/pypoetry/virtualenvs && \
    poetry --version && python3 --version && \
    poetry install --no-root --no-interaction --no-ansi

# App code last (doesn’t invalidate cfMesh)
COPY . /app
RUN find /app -type f -name "*.sh" -exec dos2unix {} \; || true && \
    poetry install --no-interaction --no-ansi

CMD ["/bin/bash"]