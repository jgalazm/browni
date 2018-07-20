import {Controller} from './Controller';
import {Model} from './Model';

let app = function(data, output, lifeCycle){

    let handler = this;
    let bathymetryReady = false;
    let initialSurfaceReady = false;
    let model, controller;


    let stringToArray = function(data){
        data = data.split('\n');
        let arr = data.map(function(row){
            return row.split(/(\s+)/).filter( function(e) {
                return e.trim().length > 0;  
            }).map(function(val){ return parseFloat(val)});; 
        });
        arr = arr.filter(function(e){
            return e.length>0; // prevents blank lines
        })
    
        return arr;
    }

    let getStringFromFile = function(url, callback){
        let xhr = new XMLHttpRequest();
        xhr.open('GET', url, true);
        xhr.responseType = 'text';

        xhr.onload = (e)=>{
            callback(xhr.responseText);
        }
        xhr.send();
    }

    let rowToMatrix = function( arr, ncols, nrows){
      let newArr = [];
      while(arr.length) newArr.push(arr.splice(0,nrows));
      return newArr;
    }

    let getArrayFromFile = function(url, callback, format='ascii'){        
        
        
      if(format == 'ascii'){
          getStringFromFile(url, (string)=>{
              let array = stringToArray(string);
              callback(array);
          });
      }
      else if(format == 'binary'){

          const xhr = new XMLHttpRequest();
          xhr.open("GET", url, true);
          xhr.responseType = 'arraybuffer';

          xhr.onload = (e)=>{
              const blob = new Blob([xhr.response], {type: "application"});
              const fileReader = new FileReader();
              let arrayBuffer;

              fileReader.onload = (event) =>{
                  arrayBuffer = event.target.result;
                  let arr = [... new Int8Array(arrayBuffer)];
                  arr = rowToMatrix( arr.slice(2,arr.length-1), arr[0], arr[1]);

                  
                  callback(arr);
  
              }
              fileReader.readAsArrayBuffer(blob);
          }

          xhr.send();
      }
    }
    
    let getArrayFromImage = function(image){
        let canvas = document.createElement('canvas');
        canvas.height = image.height;
        canvas.width = image.width;

        let ctx = canvas.getContext('2d');
        ctx.drawImage(image, 0, 0);
        let imageData = ctx.getImageData(0,0,canvas.width,canvas.height);
        
        imageData = imageData.data.filter((value, index)=>{
            return index % 4 == 0;
        });

        // convert to normal float numbers
        imageData = [...imageData];


        imageData = imageData.map((value)=>{
            return (value/255.0*(data.bathymetryMetadata.zmax-
                data.bathymetryMetadata.zmin)+
                data.bathymetryMetadata.zmin);
        });

        let matrix = [];
        while(imageData.length > 0) matrix.push(imageData.splice(0, image.width));

        return matrix;
    }

    let init = () => {
        
        this.model = new Model(data, output);

        if (lifeCycle.dataWasLoaded !== undefined){
            lifeCycle.dataWasLoaded(this.model);
        }

        this.controller = new Controller(this.model, output, lifeCycle);

        this.controller.animate();
        
    }

    let loadBathymetry = function(){
        if(data.bathymetry.slice(-3)==='png' || data.bathymetry.slice(-3)==='jpg'){
            if(!data.bathymetryMetadata){
                throw new Error('Must define data.bathymetryMetadata when using image format bathymetry');
            }
            let bathymetryImage = new Image();
            bathymetryImage.onload = function(){

                data.bathymetry = {
                    array : getArrayFromImage(bathymetryImage)
                }

                bathymetryReady = true;
                
                if(initialSurfaceReady){
                    init();
                }
            }
            bathymetryImage.src = data.bathymetry;

        }
        else{
            getArrayFromFile(data.bathymetry,function(array){
                
                data.bathymetry = {
                    array : array 
                }   
        
                bathymetryReady = true;
                
                if(initialSurfaceReady){
                    init();
                }
                
            }, data.binaryBathymetry ? 'binary':'ascii');
        }
    }

    let loadInitialCondition = function(){
        /* Detects if initialSurface, earthquake or asteroid is provided, assuming
        the user knows the right format.
        Otherwise throws an error */
        if( data.initialSurface != undefined){

            getArrayFromFile(data.initialSurface,function(array){
                
                data.initialSurface = {
                    array : array
                }
        
                initialSurfaceReady = true;
        
                if(bathymetryReady){           
                    init();
                }
                
            },'ascii');
        }
        else if( data.earthquake != undefined){
            if(typeof(data.earthquake)==="string"){
                
                getStringFromFile(data.earthquake, fileString => {
                    let earthquake = fileString.split('\n');
                    let keys = earthquake[0].split(',');
                    let key2column = {};
                    for(let i = 0; i<keys.length; i++){
                        key2column[keys[i]] = i
                    }

                    earthquake.shift();

                    earthquake = earthquake.filter((val)=> val.length>0);

                    for(let i = 0; i<earthquake.length; i++){
                        let finiteFault = earthquake[i].split(',');
                        let earthquakeDict = {};
                        keys.map((key)=>{
                            earthquakeDict[key] = (key != 'reference')? 
                                                parseFloat(finiteFault[key2column[key]]) : 
                                                finiteFault[key2column[key]]
                        });
                        earthquake[i] = earthquakeDict;
                        
                    }

                    data.earthquake = earthquake;

                    initialSurfaceReady = true;

                    if(bathymetryReady){
                        init();
                    }

                });
            }
            else{
                
                initialSurfaceReady = true;
            }

            if(bathymetryReady && initialSurfaceReady){
                init();
            }
        }
        else if(data.asteroid !== undefined){
            initialSurfaceReady = true;
            if(bathymetryReady){
                init();
            }
        }
        else {
            throw 'Need a valid data.initialSurface or data.finiteFault as input';
        }
    }

    loadBathymetry();
    loadInitialCondition();

}

export {app};
