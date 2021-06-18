FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends r-base build-essential libpq-dev libcurl4-openssl-dev python3.8 python3-setuptools python3-pip python3-dev redis
RUN R -e "install.packages(c('shinydashboard', 'reticulate', 'shiny','DT','ConfigParser','openxlsx','shinycssloaders','shinyFiles','stringr'))"
#RUN ln -s /usr/bin/python3 /usr/bin/python && \
#    ln -s /usr/bin/pip3 /usr/bin/pip
RUN apt-get update
WORKDIR /app
COPY reqs.txt reqs.txt
RUN pip install --upgrade pip
RUN pip install -r reqs.txt
RUN pip install --upgrade MODApy
RUN mkdir /usr/local/lib/python3.8/dist-packages/MODApy/logs && touch /usr/local/lib/python3.8/dist-packages/MODApy/logs/pipe_run.log
COPY ./MODApy/config.ini /usr/local/lib/python3.8/dist-packages/MODApy/
ADD start.sh .
RUN chmod +x ./start.sh
CMD ["./start.sh"]
