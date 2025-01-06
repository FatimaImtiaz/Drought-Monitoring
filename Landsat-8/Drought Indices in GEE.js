// add the study area
// var aoi: Boundary of the study area.
Map.addLayer(aoi);
Map.centerObject(aoi, 13);

// Apply cloud Mask
function maskL8sr(col) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = col.select('QA_PIXEL');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return col.updateMask(mask);
}

var startYear = 2013;
var endYear = 2021;
var startDate = ee.Date.fromYMD(startYear, 1, 1);
var endDate = ee.Date.fromYMD(endYear, 12, 31);

// Load Image collection
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
          .map(maskL8sr)
          .filterBounds(aoi)
          .filterDate(startDate, endDate);

print(l8, "l8");

// Calculate NDVI of each image
function calculateNDVI(img) {
  var ndvi = img.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI').clip(aoi);
  return ndvi.copyProperties(img, ['system:time_start', 'id']);
}

// Create NDVI image collection 
var ndvi_image_collection = l8.map(calculateNDVI);
print(ndvi_image_collection, 'NDVI Image Collection');

// Create NDVI composite for every month
var years = ee.List.sequence(startYear, endYear);
var months = ee.List.sequence(1, 12);

// Map over the years and create a monthly average collection
var monthlyImages = years.map(function(year) {
  year = ee.Number(year);  // Ensure year is an Earth Engine number
  return months.map(function(month) {
    month = ee.Number(month);  // Ensure month is an Earth Engine number
    var filtered = ndvi_image_collection
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .filter(ee.Filter.calendarRange(month, month, 'month'));
    var monthly = filtered.mean();
    return monthly.set({
      'month': month, 'year': year, 
      'system:time_start': ee.Date.fromYMD(year, month, 1).millis()
    });
  });
}).flatten();

print(monthlyImages, "Monthly NDVI Images");

// Filter out images with zero bands
var monthly_images_0band = monthlyImages.filter(ee.Filter.listContains("system:band_names", "NDVI"));
print(monthly_images_0band, 'monthly_images_0band');

// Convert the result to an ImageCollection
var monthlyCol = ee.ImageCollection.fromImages(monthlyImages);
print(monthlyCol, "Monthly NDVI Collection");

 // We now have 1 image per month for entire duratioon
var monthlyCol = ee.ImageCollection.fromImages(monthlyImages)
print (monthlyCol)


// Minimum NDVI for each month across all years
var monthlyMinImages = months.map(function(month) {
    var filtered = monthlyCol.filter(ee.Filter.eq('month', month))
    var monthlyMin = filtered.min()
    return monthlyMin.set('month', month)
})
var monthlyMin = ee.ImageCollection.fromImages(monthlyMinImages)

// Maximum NDVI for each month across all years
var monthlyMaxImages = months.map(function(month) {
    var filtered = monthlyCol.filter(ee.Filter.eq('month', month))
    var monthlyMax = filtered.max()
    return monthlyMax.set('month', month)
})
var monthlyMax = ee.ImageCollection.fromImages(monthlyMaxImages)

// Calculate VCI for 2021

//VCI for a each month in the crop season(June- September)change months and export the map.
var currentYear = 2021
var currentMonth = 6

var filtered = monthlyCol
  .filter(ee.Filter.eq('year', currentYear))
  .filter(ee.Filter.eq('month', currentMonth))
var currentMonthNdvi = ee.Image(filtered.first())

var minNdvi = ee.Image(monthlyMin.filter(ee.Filter.eq('month', currentMonth)).first())
var maxNdvi = ee.Image(monthlyMax.filter(ee.Filter.eq('month', currentMonth)).first())


var visParams = {min: 0, max: 1, palette: ['white','green']}
Map.addLayer(minNdvi, visParams, 'Minimum May NDVI')  
Map.addLayer(maxNdvi, visParams, 'Maximum May NDVI')
Map.addLayer(currentMonthNdvi, visParams, 'Current May NDVI')

// VCI = (NDVI - min) / (max - min)*100
var image = ee.Image.cat([currentMonthNdvi, minNdvi, maxNdvi]).rename(['ndvi', 'min', 'max'])
var vci = image.expression('100* (ndvi - min) / (max - min)',
    {'ndvi': image.select('ndvi'),
      'min': image.select('min'),
      'max': image.select('max')
    }).rename('vci')


var vciVisParams = {min: 0, max: 100, palette: ['white','green']}
var vciPalette = ['#a50026','#d73027','#f46d43','#fdae61',
  '#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837'];
  Map.addLayer(vci, vciVisParams, 'VCI')
  
 //Export.image.toDrive({
  image: vci,
  description: 'VCI',
  folder: 'Drought', // Specify your Google Drive folder
  region: aoi,
  scale: 30, 
  maxPixels: 1e13})


//---------------------- clculation for TCI---------------------------------
// Calculate Brightness Temperature (in Celsius)
function calculateTemp(image) {
  var temp = image.select('ST_B10').clip(aoi).multiply(0.00341802).add(149).subtract(273.15).rename('Temperature');
  return temp.copyProperties(image, ['system:time_start']);
}

var tempImageCollection = l8.map(calculateTemp);
print(tempImageCollection, 'tempImageCollection')


// Create monthly temperature composites
var years = ee.List.sequence(startYear, endYear);
var months = ee.List.sequence(1, 12);

var monthlyImages = years.map(function(year) {
  return months.map(function(month) {
    var filtered = tempImageCollection
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .filter(ee.Filter.calendarRange(month, month, 'month'));
    var monthly = filtered.mean();
    return monthly.set({'month': month, 'year': year});
  });
}).flatten();

var monthlytempCol = ee.ImageCollection.fromImages(monthlyImages);
print(monthlytempCol, 'monthlytempCol')

// Compute monthly minimum and maximum temperatures
var monthlyMinTempImages = months.map(function(month) {
  var filtered = monthlytempCol.filter(ee.Filter.eq('month', month));
  var monthlyMin = filtered.min();
  return monthlyMin.set('month', month);
});
// Minimum temperature
var monthlyMinTemp = ee.ImageCollection.fromImages(monthlyMinTempImages);

print(monthlyMinTemp, 'monthlyMinTemp')

// Maximum temperature
var monthlyMaxTempImages = months.map(function(month) {
  var filtered = monthlytempCol.filter(ee.Filter.eq('month', month));
  var monthlyMax = filtered.max();
  return monthlyMax.set('month', month);
});
var monthlyMaxTemp = ee.ImageCollection.fromImages(monthlyMaxTempImages);
print(monthlyMaxTemp, 'monthlyMaxTemp')

//TCI for a each month in the crop season(June- September)change months and export the map.
var currentYear = 2021;
var currentMonth = 6;

var filtered = monthlytempCol
  .filter(ee.Filter.eq('year', currentYear))
  .filter(ee.Filter.eq('month', currentMonth));
var currentMonthTemp = ee.Image(filtered.first());

var minTemp = ee.Image(monthlyMinTemp.filter(ee.Filter.eq('month', currentMonth)).first());
var maxTemp = ee.Image(monthlyMaxTemp.filter(ee.Filter.eq('month', currentMonth)).first());
// TCI = 100 * (LSTmax - LST) / (LSTmax – LSTmin)
var tci = currentMonthTemp.expression(
  '100 * (maxTemp - temp) / (maxTemp - minTemp)',
  {
    'temp': currentMonthTemp,
    'minTemp': minTemp,
    'maxTemp': maxTemp
  }
).rename('TCI');

var tciVisParams = {min: 0, max: 100, palette: ['blue', 'white', 'red']};
Map.addLayer(tci, tciVisParams, 'TCI');

//Export.image.toDrive({
  image: tci,
  description: 'TCI',
  folder: 'Drought', // Specify your Google Drive folder
  region: aoi,
  scale: 30, 
  maxPixels: 1e13})
//--------------------------------------------------- Calculation for  VHI---------------------------------------------------------------


var alpha = 0.5;

var vhi = vci.expression(
  'alpha * VCI + (1 - alpha) * TCI',
  {
    'VCI': vci.select('vci'),
    'TCI': tci.select('TCI'),
    'alpha': alpha
  }
).rename('VHI');

var vhiVisParams = {min: 0, max: 100, palette: ['red', 'orange', 'yellow', 'lightgreen', 'green']};
Map.addLayer(vhi, vhiVisParams, 'VHI');

//Export.image.toDrive({
  image: vhi,
  description: 'VHI',
  folder: 'Drought', // Specify your Google Drive folder
  region: aoi,
  scale: 30, 
  maxPixels: 1e13})
