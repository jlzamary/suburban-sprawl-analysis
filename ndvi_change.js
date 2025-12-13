// Vegetation change in rapidly developing suburban tech hubs? 
// A comparison of Ashburn, VA and Scottsdale, AZ (1990-2025)

// As an import in GEE (not part of the script)
var ashburn = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Point([-77.487, 39.043]);

    // Add Palettes
var palettes = require('users/gena/packages:palettes');

// Create ROI
var circle_buffer = ashburn.buffer(10000); // 10km
var roi = circle_buffer.bounds(); //400km2 total area

// Landsat 5
var ls_SR_5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');

// Landsat 9
var ls_SR_9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2");

// Cloud and scale mask
var cloudMask_LS5 = function(image) {
  var qa = image.select('QA_PIXEL');

  var mask = qa.bitwiseAnd(1 << 3).eq(0)   // shadow
             .and(qa.bitwiseAnd(1 << 4).eq(0)) // cloud
             .and(qa.bitwiseAnd(1 << 5).eq(0)); // snow

  var scaled = image.multiply(0.0000275).add(-0.2);

  return scaled.updateMask(mask)
               .copyProperties(image, image.propertyNames());
};

// Cloud mask for Landsat 9
var GetQABits = function(image, start, end, newName) {
  // set the bit extent
  var pattern = 0;
  for (var i=start; i<= end; i++) {
    pattern += Math.pow(2, i);
  }
  // Mask over all bits
  return image.select([0], [newName])
              .bitwiseAnd(pattern).rightShift(start);
};

// landsat 9 function for cloud shadow and mask
var cloudPixels = function(image){
  var QA = image.select(['QA_PIXEL']);
  return GetQABits(QA, 3, 4, 'cloud').eq(0);
};

// landsat 9 cloud confidence 
var cloud_confidence = function(image){
  var QA = image.select(['QA_PIXEL']);
  return GetQABits(QA, 10, 11, 'cloud_shadow').eq(1); // 1 = low confidence
};

// final landsat 9 mask
var maskL9clouds = function(image){
  var cp = cloudPixels(image);
  var cc = cloud_confidence(image);
  var masked = image.updateMask(cp).updateMask(cc);
  var scaled = masked.multiply(0.0000275).add(-0.2);
  return scaled.copyProperties(image, image.propertyNames());
};

var L9Unmasked = ls_SR_9.mean();
var L9masked = ls_SR_9.map(maskL9clouds).mean();


// True Color Bands (LS 9)
var ls9trueColor = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0, max: 0.3
};

// True Color Bands (LS 5)
var ls5TrueColor = {
  bands: ['SR_B3', 'SR_B2', 'SR_B1'],
  min: 0.0, max: 0.3
};

// NDVI bands (1990)
var addNDVI = function(image) {
  var ndvi = image.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI');
  return image.addBands(ndvi);
};

// NDVI bands (2025)
var addNDVI25 = function(image) {
  var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  return image.addBands(ndvi);
};


// True Color Image (1990)
var img1990 = ls_SR_5
  .filterBounds(roi)
  .filterDate('1990-06-01', '1990-09-30')
  .map(cloudMask_LS5)
  .median()
  .clip(roi);
  
// NDVI Image (1990)
var ndvi1990 = ls_SR_5
  .filterBounds(roi)
  .filterDate('1990-06-01', '1990-09-30')
  .map(cloudMask_LS5)
  .map(addNDVI)
  .median()
  .select('NDVI')
  .clip(roi);
  
// True Color Image (2025)
var img2025 = ls_SR_9
  .filterBounds(roi)
  .filterDate('2025-06-01', '2025-09-30')
  .map(maskL9clouds)
  .median()
  .clip(roi);
  
// NDVI Image (2025)
var ndvi2025 = ls_SR_9
  .filterBounds(roi)
  .filterDate('2025-06-01', '2025-09-30')
  .map(maskL9clouds)
  .map(addNDVI25)
  .median()
  .select('NDVI')
  .clip(roi);

// NDVI Color Scale 
var ndviPalette = palettes.colorbrewer.YlGn[9];
var ndviScale = {min: -1, max: 1, palette: ndviPalette};

// Change Color Scale
var changePalette = palettes.colorbrewer.PuOr[11];
var changeScale = {min: -1, max: 1, palette: changePalette};

// Calculate Change
var ndviChange = ndvi2025.subtract(ndvi1990);
print(ndviChange);

// Calculate statistics
var stats1990 = ndvi1990.reduceRegion({
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.stdDev(),
    sharedInputs: true
  }),
  geometry: roi,
  scale: 30,
  maxPixels: 1e10
});

var stats2025 = ndvi2025.reduceRegion({
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.stdDev(),
    sharedInputs: true
  }),
  geometry: roi,
  scale: 30,
  maxPixels: 1e10
});

print('NDVI Statistics');
print('1990 Mean NDVI:', stats1990.get('NDVI_mean'));
print('1990 StdDev:', stats1990.get('NDVI_stdDev'));
print('2025 Mean NDVI:', stats2025.get('NDVI_mean'));
print('2025 StdDev:', stats2025.get('NDVI_stdDev'));


// Calculate vegetation loss and gain
var threshold = 0.1;

// Loss
var loss = ndviChange.lt(-threshold).selfMask();
var lossArea = loss.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: roi,
  scale: 30,
  maxPixels: 1e10
});

// Gain
var gain = ndviChange.gt(threshold).selfMask();
var gainArea = gain.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: roi,
  scale: 30,
  maxPixels: 1e10
});

print('Total Area Change Statistics');
print('Vegetation Loss (km2):', ee.Number(lossArea.get('NDVI')).divide(1e6));
print('Vegetation Gain (km2):', ee.Number(gainArea.get('NDVI')).divide(1e6));

// Net change
var netChange = ee.Number(gainArea.get('NDVI')).subtract(ee.Number(lossArea.get('NDVI'))).divide(1e6);
print('Net Vegetation Change (km2):', netChange);

// Add map layers
Map.addLayer(img1990, ls5TrueColor, "1985 True Color");
Map.addLayer(ndvi1990, {min: -0.2, max: 1, palette: ['white', 'black']}, "Gray Scale 90");
Map.addLayer(ndvi1990, ndviScale, "NDVI 1990");

Map.addLayer(img2025, ls9trueColor, '2025 True Color');
Map.addLayer(ndvi2025, {min: -0.2, max: 1, palette: ['white', 'black']}, "Gray Scale 25");
Map.addLayer(ndvi2025, ndviScale, "NDVI 2025");

Map.addLayer(ndviChange, changeScale, "NDVI Change");

Map.centerObject(roi);