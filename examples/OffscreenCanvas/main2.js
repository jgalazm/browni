let canvas = document.getElementById('canvas');
var offscreen = canvas.transferControlToOffscreen();

var worker = new Worker("offscreen.js"); 
worker.postMessage({canvas: offscreen}, [offscreen]);