FROM node:16.0.0-alpine3.13 as builder

# Install and cache app dependencies
WORKDIR /nami
COPY package.json .
COPY package-lock.json .
RUN yarn install

COPY babel.config.js .
COPY rollup.config.js .
COPY src src

RUN yarn build

# copy compiled file to smaller image
FROM alpine:3.10
COPY --from=builder /nami/build /nami/build

VOLUME /nami/build