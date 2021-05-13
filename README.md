# Nami

WebGL tsunami simulator that seamlessly integrates with the web.

See [this demo](https://codepen.io/jgalazm/pen/MRaWVL): 

# Installation 

Nami can be installed using the Node Package Manager

```
npm install namijs
```

* Import with <script/> (coming soon)
* Import with 'import' (coming soon)
* Import with 'require' (coming soon)

# Examples 
* ["Hello world"](https://codepen.io/jgalazm/pen/MRaWVL)
* Switch scenario with key press (coming soon)
* Integration with Three.js, Leaflet, ArcGIS and CesiumJS (coming soon)
* Save colored texture to a video (coming soon)

# API
(coming soon)

# More info
- Why Nami?
- How was it done? 
- Has it been used before? 
- Are simulations accurate? 
- Can I simulate "this"?

Check [the preprint](hal.inria.fr/hal-02112763)

# Local development environment

A 'docker-compose.yml' is provided with two services `nami` and `local-server`. The `nami` service opens a container whose image, once built, contains an updated-from-source build of the library. This build is created with rollup and stored in a docker volume called `build-volume` which is shared with the `local-server`. This `local-server` is a python simple http server that will enable properly loading this built file.