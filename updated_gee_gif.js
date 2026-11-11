// Day 09 – Urban/Rural: Ashburn, VA summer NDVI, 1990–2024 (water masked)
// YouTube Video: https://www.youtube.com/watch?v=Wsfis3IxBtc
// Original Code: https://code.earthengine.google.com/d1c85a213097e0ba4b6f53fdf9941181
//
// Center point: Ashburn, VA (Loudoun County)
//
// AN UPDATE OF THE ORIGINAL SCRIPT WITH THE HELP OF CLAUDE CODE
//

// ===================== Settings =====================
var BUFFER_M        = 10000;  // 10 km radius -> ~20 x 20 km box
var YEARS           = [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024];
var SEASON_START    = '-06-01';
var SEASON_END      = '-09-30';

var MNDWI_THRESHOLD = 0;      // composite MNDWI above this = water (raise to ~0.1 if rooftops or bare fields get flagged)
var WATER_BUFFER_M  = 0;     // grow water mask to drop mixed shoreline / stream-edge pixels (0 = off)
var HARMONIZE_OLI   = true;   // adjust Landsat 8/9 reflectance to TM/ETM+ (Roy et al. 2016)

// ===================== ROI =====================
var ashburn = ee.Geometry.Point([-77.487, 39.043]);
var roi = ashburn.buffer(BUFFER_M).bounds();

// ===================== Collections =====================
var ls5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
var ls8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var ls9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');

// Common band names so every sensor goes through the same code
var BANDS_OUT = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
var TM_BANDS  = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']; // Landsat 5
var OLI_BANDS = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']; // Landsat 8/9

// OLI -> ETM+ surface reflectance harmonization (Roy et al. 2016)
var OLI_SLOPE = ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]);
var OLI_ITCP  = ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]);

// ===================== Masking & scaling =====================
// QA_PIXEL bits: 0 fill, 1 dilated cloud, 2 cirrus, 3 cloud, 4 cloud shadow, 5 snow
function maskQA(image) {
  var qa = image.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(parseInt('111111', 2)).eq(0);
  return image.updateMask(clear);
}

// Landsat 5: mask, keep optical bands only, scale to reflectance
function prepTM(image) {
  var sr = maskQA(image)
    .select(TM_BANDS, BANDS_OUT)
    .multiply(0.0000275).add(-0.2);
  return ee.Image(sr.copyProperties(image, ['system:time_start']));
}

// Landsat 8/9: same, plus optional harmonization to TM/ETM+
function prepOLI(image) {
  var sr = maskQA(image)
    .select(OLI_BANDS, BANDS_OUT)
    .multiply(0.0000275).add(-0.2);
  if (HARMONIZE_OLI) {
    sr = sr.multiply(OLI_SLOPE).add(OLI_ITCP);
  }
  return ee.Image(sr.copyProperties(image, ['system:time_start']));
}

// ===================== Yearly composites =====================
function seasonalComposite(year) {
  var start = year + SEASON_START;
  var end   = year + SEASON_END;

  var col;
  if (year <= 2011) {
    col = ls5.filterBounds(roi).filterDate(start, end).map(prepTM);
  } else {
    // Landsat 8 + 9 together (Landsat 9 only adds scenes from late 2021 on)
    col = ls8.merge(ls9).filterBounds(roi).filterDate(start, end).map(prepOLI);
  }

  var composite = col.median();

  var ndvi  = composite.normalizedDifference(['nir', 'red']).rename('NDVI');
  var mndwi = composite.normalizedDifference(['green', 'swir1']).rename('MNDWI');

  // One water mask per year, built from the composite
  var water = mndwi.gt(MNDWI_THRESHOLD);
  if (WATER_BUFFER_M > 0) {
    water = water.focalMax({radius: WATER_BUFFER_M, units: 'meters', kernelType: 'square'});
  }
  water = water.rename('water');

  return ndvi.updateMask(water.not())   // NDVI over land only
    .addBands(water)
    .clip(roi)
    .set('year', year)
    .set('n_scenes', col.size());
}

var ndviCollection = ee.ImageCollection.fromImages(YEARS.map(seasonalComposite));

// ===================== Visualization =====================
var ndviVis = {
  min: 0,
  max: 0.85,
  palette: ['f2efe6', 'd9c89e', 'a8c66c', '4a8c2a', '0b4d12'] // bare/urban -> lawns/fields -> forest
};

var noDataRGB = ee.Image.constant([158, 158, 158]); // gray = no clear observations
var waterRGB  = ee.Image.constant([47, 111, 191]);  // blue = masked water

var frames = ndviCollection.map(function(img) {
  var background = noDataRGB.visualize({min: 0, max: 255});
  var water = waterRGB.updateMask(img.select('water').selfMask()).visualize({min: 0, max: 255});
  var veg = img.select('NDVI').visualize(ndviVis);
  return background.blend(water).blend(veg).clip(roi);
});

var gifParams = {
  region: roi,
  dimensions: 600,
  crs: 'EPSG:3857',
  framesPerSecond: 1
};

// ===================== Display =====================
print('NDVI Animation GIF URL:');
print(frames.getVideoThumbURL(gifParams));
print(ui.Thumbnail(frames, gifParams));

// Check the water mask on the latest year (Inspector: Potomac, Goose Creek, Beaverdam Reservoir, ponds)
var latest = ndviCollection.filter(ee.Filter.eq('year', 2024)).first();
Map.centerObject(roi, 12);
Map.addLayer(latest.select('NDVI'), ndviVis, 'NDVI 2024 (land only)');
Map.addLayer(latest.select('water').selfMask(), {palette: ['2f6fbf']}, 'Water mask 2024');

// Quick stats: land-only mean NDVI and scene count per year
var stats = ee.FeatureCollection(ndviCollection.map(function(img) {
  var meanNdvi = img.select('NDVI').reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 30,
    maxPixels: 1e13
  }).get('NDVI');
  return ee.Feature(null, {
    year: img.get('year'),
    mean_land_NDVI: meanNdvi,
    n_scenes: img.get('n_scenes')
  });
}));
print('Per-year stats', stats);
print(ui.Chart.feature.byFeature(stats, 'year', ['mean_land_NDVI'])
  .setOptions({title: 'Mean summer NDVI (water excluded)', hAxis: {format: '####'}}));

// ===================== Export =====================
Export.video.toDrive({
  collection: frames,
  description: 'Ashburn_NDVI_GIF',
  folder: 'EarthEngine',
  fileNamePrefix: 'Ashburn_NDVI_GIF',
  region: roi,
  scale: 30,
  crs: 'EPSG:3857',
  framesPerSecond: 1,
  maxPixels: 1e13
});