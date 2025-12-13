// YouTube Video: https://www.youtube.com/watch?v=Wsfis3IxBtc
// Original Code: https://code.earthengine.google.com/d1c85a213097e0ba4b6f53fdf9941181

// As an import in GEE (not part of the script)
var ashburn = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Point([-77.487, 39.043]);

// Create ROI
var circle_buffer = ashburn.buffer(10000); // 10km
var roi = circle_buffer.bounds(); //400km2 total area

// Landsat Collection
var ls5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
var ls8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var ls9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');

// Scale
function scaleSR(image) {
  return image.multiply(0.0000275).add(-0.2)
              .copyProperties(image, image.propertyNames());
}

// QA Bit 
function getQABits(image, start, end) {
  var pattern = 0;
  for (var i = start; i <= end; i++) {
    pattern += Math.pow(2, i);
  }
  return image
    .bitwiseAnd(pattern)
    .rightShift(start);
}

// Cloud, Shadow, and Snow mask 
function maskCloudsC2(image) {
  var qa = image.select('QA_PIXEL');

  var cloudShadow = getQABits(qa, 3, 3).eq(0);
  var cloud       = getQABits(qa, 4, 4).eq(0);
  var snow        = getQABits(qa, 5, 5).eq(0);

  return image
    .updateMask(cloudShadow)
    .updateMask(cloud)
    .updateMask(snow);
}

// Calculate NDVI

// Landsat 5
function addNDVI_L5(image) {
  return image.addBands(
    image.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI')
  );
}

// Landsat 8/9
function addNDVI_L89(image) {
  return image.addBands(
    image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')
  );
}

// Mask water with NDWI
// Landsat 5
function addNDWI_L5(image) {
  return image.addBands(
    image.normalizedDifference(['SR_B2', 'SR_B4']).rename('NDWI')
  );
}

// Landsat 8/9
function addNDWI_L89(image) {
  return image.addBands(
    image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI')
  );
}

function maskWater(image) {
  return image.updateMask(image.select('NDWI').lt(0.2));
}

// Select years
var years = [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024];

// Images per year
var ndviCollection = ee.ImageCollection.fromImages(
  years.map(function(year) {
    var start = year + '-06-01';
    var end   = year + '-09-30';
    var image;

    if (year <= 2011) {
      image = ls5
        .filterBounds(roi)
        .filterDate(start, end)
        .map(maskCloudsC2)
        .map(scaleSR)
        .map(addNDVI_L5)
        .map(addNDWI_L5)
        .map(maskWater)
        .median();

    } else if (year <= 2021) {
      image = ls8
        .filterBounds(roi)
        .filterDate(start, end)
        .map(maskCloudsC2)
        .map(scaleSR)
        .map(addNDVI_L89)
        .map(addNDWI_L89)
        .map(maskWater)
        .median();

    } else {
      image = ls9
        .filterBounds(roi)
        .filterDate(start, end)
        .map(maskCloudsC2)
        .map(scaleSR)
        .map(addNDVI_L89)
        .map(addNDWI_L89)
        .map(maskWater)
        .median();
    }

    return image
      .select('NDVI')
      .clip(roi)
      .set('year', year);
  })
);

// Choose palettw
var ndviVis = {
  min: -0.2,
  max: 0.8,
  palette: ['white', 'yellow', 'green', 'darkgreen']
};

// GIF frames and settings
var frames = ndviCollection.map(function(img) {
  return img.visualize(ndviVis).clip(roi);
});

var gifParams = {
  region: roi,
  dimensions: 600,
  crs: 'EPSG:3857',
  framesPerSecond: 1
};

// Display
print('NDVI Animation GIF URL:');
print(frames.getVideoThumbURL(gifParams));
print(ui.Thumbnail(frames, gifParams));

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
