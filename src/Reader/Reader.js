import {getArrayFromImage} from './FileUtils';

const Reader = function(data, init){

    let loadBathymetry = function(resolve){

        if(typeof data.bathymetry == 'object'){
            // assume bathymetry is an array
            data.bathymetry = {
                array: data.bathymetry
            };
            resolve(data.bathymetry);
        }
        else if(data.bathymetry.slice(-3)==='png' || data.bathymetry.slice(-3)==='jpg'){
            if(!data.bathymetryMetadata){
                throw new Error('Must define data.bathymetryMetadata when using image format bathymetry');
            }
            let bathymetryImage = new Image();
            bathymetryImage.onload = () =>{

                data.bathymetry = {
                    array : getArrayFromImage(bathymetryImage, data.bathymetryMetadata)
                }

                resolve(data.bathymetry.array);
            }
            bathymetryImage.src = data.bathymetry;

        }
        else{
            getArrayFromFile(data.bathymetry, (array) =>{
                
                data.bathymetry = {
                    array : array 
                }   
        
                resolve(data.bathymetry.array);
                
            }, data.binaryBathymetry ? 'binary':'ascii');
        }
    }

    let loadInitialCondition = function(resolve){
        /* Detects if initialSurface, earthquake or asteroid is provided, assuming
        the user knows the right format.
        Otherwise throws an error */

        if( data.initialSurface != undefined){

            getArrayFromFile(data.initialSurface.file,function(array){
                
                data.initialSurface.array = array;
        
                initialSurfaceReady = true;
                
                resolve(data.initialSurface.array);
                
            },'ascii');
        }
        else if( data.earthquake != undefined){
            if(typeof(data.earthquake)==="string"){
                
                getStringFromFile(data.earthquake, fileString => {
                    let earthquake = fileString.replace(String.fromCharCode(13),'').split('\n');
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

                    resolve(earthquake);
                });
            }
            else{
                
                resolve(data.earthquake);
            }

            if(bathymetryReady && initialSurfaceReady){
                init();
            }
        }
        else if(data.asteroid !== undefined){
            if(bathymetryReady){
                init();
            }
        }
        else {
            throw 'Need a valid data.initialSurface or data.finiteFault as input';
        }
    }

  
        

    let newData = {};

    newData.bathymetry = new Promise( (resolve, reject) =>{
        loadBathymetry(resolve);
    });

    newData.initialCondition = new Promise( (resolve, reject) => {
        loadInitialCondition(resolve);
    });
    
    return newData;
}

export default Reader;