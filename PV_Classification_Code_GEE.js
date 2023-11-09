/**
* Function to mask clouds using the Sentinel-2 QA band
* @param {ee.Image} image Sentinel-2 image
* @return {ee.Image} cloud masked Sentinel-2 image
*/
function maskS2clouds(image) {
  var qa = image.select('QA60'); 
  
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
      
  return image.updateMask(mask).divide(10000);
}

var roi = table.geometry();
Map.addLayer(roi);
var image = ee.ImageCollection('COPERNICUS/S2_SR')
              .filterBounds(roi)
              .filterDate('2020-02-01', '2020-05-01')
              // Pre-filter to get less cloudy granules.
              .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
              .map(maskS2clouds)
              .select('B8','B4','B2','B3','B11','B12').median();
var rgbVis = {
  min: 0.0,
  max: 0.3,
  gamma: 1.4,
  bands: ['B4', 'B3', 'B2'],
};
Map.addLayer(image.clip(roi), rgbVis, 'RGB_0807');

/// Elevation
var DEMdataset = ee.Image('JAXA/ALOS/AW3D30_V1_1').clip(roi);
var elevation = DEMdataset.select('AVE');
// Slope
var slope = ee.Terrain.slope(elevation).clip(roi);
// Aspect
var aspect = ee.Terrain.aspect(elevation);
var sinImage = aspect.divide(180).multiply(Math.PI).sin();
//Spectral Index
function cal_NDVI(data){
  var ndvi0 = data.expression('(nir-red)/(nir+red)',{'nir' :data.select('B8'),'red' :data.select('B4')});
  var ndvi = ee.Image(ndvi0);
  return ndvi;
}
function cal_MNDWI(data){
  var ndwi0 = data.expression('(green-mir)/(mir+green)',{'mir' :data.select('B11'),'green' :data.select('B3')});
  var mndwi = ee.Image(ndwi0);
  return mndwi;
}
function cal_NDPI(data){ 
  var ndpi0 = data.expression('(swir1-nir)/(nir-swir2)',{'nir' :data.select('B8'),'swir1' :data.select('B11'),'swir2' :data.select('B12')});
  var ndpi = ee.Image(ndpi0);
  return ndpi;}
function cal_NDBI(data){
  var ndbi0 = data.expression('(swir1-nir)/(nir+swir1)',{'nir' :data.select('B8'),'swir1' : data.select('B11')});
  var ndbi = ee.Image(ndbi0);
  return ndbi;
}
function cal_SAVI(data){
  var savi0 = data.expression('1.5*(nir-red)/(nir+red+0.5)',{'nir' :data.select('B8'),'red' :data.select('B4')});
  var savi = ee.Image(savi0);
  return savi;
}

var ndvi = cal_NDVI(image).rename('ndvi');
var mndwi = cal_MNDWI(image).rename('mndwi');
var ndpi = cal_NDPI(image).rename('ndpi');
var ndbi = cal_NDBI(image).rename('ndbi');
var savi = cal_SAVI(image).rename('savi');

//GLCM
var getGLCM = function(data) {
  var asm = data.int32().glcmTexture({size: 4}).select('B8_asm');
  var contrast = data.int32().glcmTexture({size: 4}).select('B8_contrast');
  var corr = data.int32().glcmTexture({size: 4}).select('B8_corr');
  var vari = data.int32().glcmTexture({size: 4}).select('B8_var');
  var idm = data.int32().glcmTexture({size: 4}).select('B8_idm');
  var savg = data.int32().glcmTexture({size: 4}).select('B8_savg');
  var ent = data.int32().glcmTexture({size: 4}).select('B8_ent');
  var diss = data.int32().glcmTexture({size: 4}).select('B8_diss');
  return data.addBands(asm).addBands(contrast).addBands(corr).addBands(vari).addBands(idm).addBands(ent).addBands(savg).addBands(diss);
};

var glcm = getGLCM(image).select('B8_asm','B8_contrast','B8_corr','B8_var','B8_idm','B8_savg','B8_ent','B8_diss');

var dataset = image.addBands(ndvi).addBands(mndwi).addBands(ndpi).addBands(ndbi).addBands(savi).addBands(glcm)
              .addBands(elevation).addBands(slope).addBands(sinImage);

// Create 40000 random points in the region.
var randomPointsphotovoltaic = ee.FeatureCollection.randomPoints(photovoltaic,40000,0);
var randomPointsother = ee.FeatureCollection.randomPoints(other,40000,0);

var addproperty0 = function(feature) {
  return feature.set('landcover',0);
};
var addproperty1 = function(feature) {
  return feature.set('landcover',1);
};

// Map the area getting function over the FeatureCollection.
var pointother  = randomPointsother.map(addproperty0);
var pointphotovoltaic = randomPointsphotovoltaic.map(addproperty1);
var pointother = pointother.randomColumn();
var pointphotovoltaic = pointphotovoltaic.randomColumn();

// Split Training data and Validation data
var split = 0.7;  // 7:3
// photovoltaic
var Training_data_photovoltaic = pointphotovoltaic.filter(ee.Filter.lte('random', split));
print('Training_data_photovoltaic',Training_data_photovoltaic)
var Validation_data_photovoltaic = pointphotovoltaic.filter(ee.Filter.gte('random', split));
print('Validation_data_photovoltaic',Validation_data_photovoltaic)
// other
var Training_data_other = pointother.filter(ee.Filter.lte('random', split));
print('Training_data_other',Training_data_other)
var Validation_data_other = pointother.filter(ee.Filter.gte('random', split));
print('Validation_data_other',Validation_data_other)
//merging points
var Training_points = Training_data_photovoltaic.merge(Training_data_other);
var Validation_points = Validation_data_photovoltaic.merge(Validation_data_other);

var training = dataset.sampleRegions({
  collection: Training_points,
  properties: ['landcover'], 
  scale:30,
  tileScale: 16
});

// Make a Random Forest classifier and train it.
var classifier = ee.Classifier.smileRandomForest(80,7).train(training,'landcover');

// Classify the image.
var classified = dataset.clip(roi).classify(classifier);
var kernel = ee.Kernel.square({radius: 1});
// Perform an erosion followed by a dilation, display.
var opened = classified
            .focal_min({kernel: kernel, iterations: 1})
            .focal_max({kernel: kernel, iterations: 1});

var palette = [
'#000000',    // 'photovoltaic(1)'
'#ffffff',    // 'other(0)'
];


// Map.addLayer(classified, {min: 0, max: 1, palette: palette}, 'photovoltaic Classification');
// Map.addLayer(opened, {min: 0, max: 1, palette: palette}, 'DP OPENED');


//Sample the input imagery to get a FeatureCollection of Validation data.
var Validation_points = dataset.sampleRegions({
  collection: Validation_points,
  properties: ['landcover'],
  scale: 30,
  tileScale: 16
});

// Classify the validation data.
var validated = Validation_points.classify(classifier,'testlandcover');
 
// Get a confusion matrix representing expected accuracy.
var testAccuracy = validated.errorMatrix('landcover','testlandcover');
//Printing of confusion matrix may time out. 
//Alternatively, you can export it as CSV
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());

//Calculate and print the following assessment metrics. 
//Producer's accuracy, Consumer's accuracy, Kappa coefficient
var producersAccuracy = testAccuracy.producersAccuracy();
var consumerAccuracy = testAccuracy.consumersAccuracy();
var kappaCoefficient = testAccuracy.kappa();
print('producersAccuracy:', producersAccuracy);
print('consumerAccuracy:', consumerAccuracy);
print('kappaCoefficient:', kappaCoefficient);

  
