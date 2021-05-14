FROM node:16.0.0-alpine3.13 as builder

# Install and cache app dependencies
WORKDIR /nami
COPY package.json .
COPY package-lock.json .
RUN yarn install

CMD yarn build