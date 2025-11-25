FROM ubuntu:22.04

USER root

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && \
    apt-get install -y python3.11 python3.11-venv python3-pip gnupg wget curl software-properties-common python3-apt

RUN wget -q -O - https://dl.openfoam.org/gpg.key | apt-key add -
RUN add-apt-repository http://dl.openfoam.org/ubuntu
RUN apt-get update
RUN apt-get install -y openfoam11

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install poetry

WORKDIR /app

COPY . /app

RUN poetry install --no-interaction --no-ansi

CMD ["/bin/bash"]