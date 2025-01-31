FROM ubuntu
FROM python:3.7-alpine

ENV ZHUNT_HOME=/zhunt
ARG REQUIREMENTS=requirements.txt

RUN mkdir -p $ZHUNT_HOME/bin\
    mkdir -p $ZHUNT_HOME/test\
    mkdir -p $ZHUNT_HOME/src

ENV FLASK_APP app.py
ENV FLASK_RUN_HOST 0.0.0.0
RUN apk add --no-cache gcc make musl-dev

# python dependency installs
COPY $REQUIREMENTS $ZHUNT_HOME/$REQUIREMENTS
RUN pip install -r $ZHUNT_HOME/$REQUIREMENTS

# path changes
ENV PYTHONPATH $PYTHONPATH:$ZHUNT_HOME
ENV PATH $PATH:$ZHUNT_HOME/bin

# add source code and tests
COPY test $ZHUNT_HOME/test
COPY src $ZHUNT_HOME/src

# compile zhunt
RUN make -C $ZHUNT_HOME/src TARGET=$ZHUNT_HOME/bin/zhunt

WORKDIR $ZHUNT_HOME

COPY . .
CMD ["flask", "run"]
