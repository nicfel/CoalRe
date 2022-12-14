# Dockerfile to build container for unit testing

FROM openjdk:8

RUN apt-get update && apt-get install -y git ant

WORKDIR /root

ADD . ./

ENTRYPOINT ant test



# Dockerfile to build container for unit testing

FROM debian:stable

RUN apt-get update
RUN apt-get install -y openjdk-17-jdk openjfx
RUN apt-get install -y ant
RUN apt-get install -y jblas

WORKDIR /root

ADD . ./


ENTRYPOINT JAVA_FX_HOME=/usr/share/java/ ant test