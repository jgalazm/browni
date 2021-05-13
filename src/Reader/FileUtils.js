export function stringToArray(data) {
  data = data.split("\n");
  let arr = data.map(function(row) {
    return row
      .split(/(\s+)/)
      .filter(function(e) {
        return e.trim().length > 0;
      })
      .map(function(val) {
        return parseFloat(val);
      });
  });
  arr = arr.filter(function(e) {
    return e.length > 0; // prevents blank lines
  });

  return arr;
}

export function getStringFromFile(url, callback) {
  let xhr = new XMLHttpRequest();
  xhr.open("GET", url, true);
  xhr.responseType = "text";

  xhr.onload = e => {
    callback(xhr.responseText);
  };
  xhr.send();
}

export function rowToMatrix(arr, ncols, nrows) {
  let newArr = [];
  while (arr.length) newArr.push(arr.splice(0, nrows));
  return newArr;
}

export function getArrayFromFile(url, callback, format = "ascii") {
  if (format == "ascii") {
    getStringFromFile(url, string => {
      let array = stringToArray(string);
      callback(array);
    });
  } else if (format == "binary") {
    const xhr = new XMLHttpRequest();
    xhr.open("GET", url, true);
    xhr.responseType = "arraybuffer";

    xhr.onload = e => {
      const blob = new Blob([xhr.response], { type: "application" });
      const fileReader = new FileReader();
      let arrayBuffer;

      fileReader.onload = event => {
        arrayBuffer = event.target.result;
        let arr = [...new Float64Array(arrayBuffer)];

        callback(arr);
      };
      fileReader.readAsArrayBuffer(blob);
    };

    xhr.send();
  }
}

export function getArrayFromImage(image, bathymetryMetadata) {
  image.crossOrigin = "Anonymous";
  
  let canvas = document.createElement("canvas");
  canvas.height = image.height;
  canvas.width = image.width;

  let ctx = canvas.getContext("2d");
  ctx.drawImage(image, 0, 0);
  let imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);

  imageData = imageData.data.filter((value, index) => {
    return index % 4 == 0;
  });

  // convert to normal float numbers
  imageData = [...imageData];

  imageData = imageData.map(value => {
    return (
      (value / 255.0) * (bathymetryMetadata.zmax - bathymetryMetadata.zmin) +
      bathymetryMetadata.zmin
    );
  });

  let matrix = [];
  while (imageData.length > 0) matrix.push(imageData.splice(0, image.width));

  return matrix;
}
