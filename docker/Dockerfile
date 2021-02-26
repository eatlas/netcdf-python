#
# Docker image useful for development within this project.
#

# Based on latest Python 3.
FROM python:3

# Create a new user and a data directory, and make the new user the owner of the data directory.
RUN groupadd -r python && useradd -r -u 1000 -g python python \
    && mkdir /data && chown -R python:python /data

# Update the OS and install any utilities that are required.
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq \
    && apt-get upgrade -y \
    && apt-get install libproj-dev libgeos-dev -y \
    && pip3 install --upgrade pip numpy

# Install dependent Python packages.
COPY requirements.txt .
RUN pip3 --no-cache-dir install -r requirements.txt

# Use the new user.
USER python
