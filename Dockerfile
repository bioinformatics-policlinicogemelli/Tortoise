FROM continuumio/miniconda3:25.1.1-2
#IMPORT FILES
COPY assets             /tortoise/assets
COPY lib                /tortoise/lib
COPY tortoise.py        /tortoise
COPY main.py            /tortoise
COPY environment.yaml   /tortoise
#UPDATE AND INSTALL APT PACKAGES
RUN apt update && apt install build-essential -y
RUN conda env create -n tortoise --file=/tortoise/environment.yaml
#START
WORKDIR /tortoise
EXPOSE 8593
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "tortoise", "python", "main.py"]
