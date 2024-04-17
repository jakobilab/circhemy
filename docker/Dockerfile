FROM ubuntu:22.04

LABEL stage=builder

LABEL maintainer="tjakobi@arizona.edu"

ARG CIRCHEMY_VERSION

ARG BUILD_DATE

LABEL org.label-schema.build-date=$BUILD_DATE

# update all repos, install minimal start packages

RUN sed -i 's/http:\/\/archive\.ubuntu\.com/http:\/\/mirror\.arizona\.edu/' /etc/apt/sources.list; sed -i 's/http:\/\/security\.ubuntu\.com/http:\/\/mirror\.arizona\.edu/' /etc/apt/sources.list; apt-get update; apt-get install --no-install-recommends wget git gpg ca-certificates -y

RUN useradd -ms /bin/bash circhemy

# install R and Python

RUN apt-get update; DEBIAN_FRONTEND=noninteractive TZ=America/Phoenix apt-get install --no-install-recommends python3 python3-dev bzip2 -y;

# clone circhemy
# install pip
# install circhemy
# unpack database (as the container is immutable)

RUN git config --global url."https://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@github.com/".insteadOf "https://github.com/"

RUN git clone https://github.com/jakobilab/circhemy.git ; ls -la circhemy/;

RUN wget https://bootstrap.pypa.io/get-pip.py; python3 get-pip.py; rm get-pip.py; python3 -m pip install -U setuptools; python3 -m pip install circhemy/; pip cache purge; apt-get purge python3-dev; apt-get autoremove; apt-get autoclean; rm -rf /var/lib/apt/lists/*; rm circhemy/ -rf; bunzip2 /usr/local/lib/python3.10/dist-packages/circhemy/data/circhemy.sqlite3.bz2;

USER circhemy

# define entrypoint

ENTRYPOINT ["circhemy"]
