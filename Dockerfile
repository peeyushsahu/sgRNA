FROM python:3.8

# Set some meta-data about this container
LABEL maintainer="Peeyush Sahu"

# install java and update
RUN apt update && \
    apt install -y openjdk-11-jdk && \
    apt-get autoremove -y && \
    apt-get clean

# Create work directory
RUN mkdir -p /app/tools
WORKDIR /app

# Install nextflow
RUN wget -qO- https://get.nextflow.io | bash
RUN chmod +x nextflow

# Create other required folders and upload needed files
RUN mkdir /app/tools
COPY tools/bowtie2-2.5.1-linux-x86_64.zip /app/tools/
RUN unzip /app/tools/bowtie2-2.5.1-linux-x86_64.zip -d /app/tools/
RUN rm /app/tools/bowtie2-2.5.1-linux-x86_64.zip

COPY nextflow.config /app/
COPY main.nf /app/
COPY sgRNA_alignment.nf /app/
COPY sgRNA_processing.py /app/
RUN mkdir /app/data

# Install python package
RUN pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org --upgrade pip
RUN pip install pandas==2.0.3
#RUN pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org --no-cache-dir -r requirements.txt

# Set entry point
ENTRYPOINT ./nextflow -c nextflow.config run main.nf