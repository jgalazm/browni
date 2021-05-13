FROM node:16.0.0-alpine3.13 as builder

# Install and cache app dependencies
WORKDIR /nami
COPY package.json .
COPY package-lock.json .
RUN yarn install

COPY babel.config.js .
COPY rollup.config.js .
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
COPY src src

RUN yarn build

# copy compiled file to smaller image
# FROM alpine:3.10
# COPY --from=builder /nami/build /nami/build

VOLUME /nami/build