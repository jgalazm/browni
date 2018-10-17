var map = L.map('map').setView([37.8, -96], 4);

L.tileLayer('https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token={accessToken}', {
    attribution: 'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, <a href="https://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, Imagery © <a href="https://www.mapbox.com/">Mapbox</a>',
    maxZoom: 18,
    id: 'mapbox.streets',
    accessToken: 'pk.eyJ1Ijoiam9zZWdhbGF6IiwiYSI6ImNqbmRrd3cyajA0YTUzcHBnMWNrM3M5ZXkifQ.ZiavFKiPD__VnlEG3QLrFw'
}).addTo(map);

var bounds = L.latLngBounds([[32, -130], [13, -100]]);

var videoElement = document.getElementById('videoTarget');

var videoOverlay = L.videoOverlay(videoElement, bounds, {
    opacity: 1.0
}).addTo(map);   