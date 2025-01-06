
// Load image collections and study area boundary
var modisLST:Imagecollection MODIS/061/MOD11A2
var modis:Imagecollection MODIS/061/MOD13Q1
var aoi: Boundary of the study area.

var startYear = 2012;
var endYear = 2021;
var startDate = ee.Date.fromYMD(startYear, 1, 1);
var endDate = ee.Date.fromYMD(endYear, 12, 31);
var filtered = modis
  .filter(ee.Filter.date(startDate, endDate))
  .filter(ee.Filter.bounds(aoi));

// Cloud Masking
var bitwiseExtract = function(input, fromBit, toBit) {
  var maskSize = ee.Number(1).add(toBit).subtract(fromBit);
  var mask = ee.Number(1).leftShift(maskSize).subtract(1);
  return input.rightShift(fromBit).bitwiseAnd(mask);
};

var maskSnowAndClouds = function(image) {
  var summaryQa = image.select('SummaryQA');
  var qaMask = bitwiseExtract(summaryQa, 0, 1).lte(1);
  var maskedImage = image.updateMask(qaMask);
  return maskedImage.copyProperties(image, ['system:time_start']);
};

var maskedCol = filtered.map(maskSnowAndClouds);

// selection of NDVI band
var ndviCol = maskedCol.select('NDVI');
var scaleNdvi = function(image) {
  var scaled = image.divide(10000);
  return scaled.copyProperties(image, ['system:time_start']);
};
// Scaling
var ndviScaled = ndviCol.map(scaleNdvi);

// Create NDVI composite for every month
var createMonthlyComposite = function(date) {
  var year = ee.Date(date).get('year');
  var month = ee.Date(date).get('month');
  var filtered = ndviScaled
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(month, month, 'month'));
  var monthly = filtered.mean();
  return monthly.set({
    'system:time_start': ee.Date(date).millis(),
    'month': month,
    'year': year
  });
};

var monthlyDates = ee.List.sequence(0, (endYear - startYear + 1) * 12 - 1)
  .map(function(monthOffset) {
    return startDate.advance(monthOffset, 'month');
  });

var monthlyCol = ee.ImageCollection(monthlyDates.map(createMonthlyComposite));
print(monthlyCol, 'monthly');

// Compute Minimum and Maximum NDVI for each month across all years
var monthlyMinMax = ee.List.sequence(1, 12).map(function(month) {
  var filtered = monthlyCol.filter(ee.Filter.eq('month', month));
  var min = filtered.min().set('month', month);
  var max = filtered.max().set('month', month);
  return ee.Feature(null, {month: month, min: min, max: max});
});
//---------------------------------------------------- Calculation of VCI ----------------------------------------
// Calculate VCI for each image in the monthly collection
var calculateVCI = function(image) {
  var month = ee.Number(image.get('month'));
  var minMax = ee.Feature(monthlyMinMax.filter(ee.Filter.eq('month', month)).get(0));
  var min = ee.Image(minMax.get('min'));
  var max = ee.Image(minMax.get('max'));
  var vci = image.subtract(min).divide(max.subtract(min)).multiply(100).rename('VCI'); // Multiply by 100 for scaling
  return image.addBands(vci).select('VCI').copyProperties(image, ['system:time_start', 'month', 'year']);
};

var vciCol = monthlyCol.map(calculateVCI);
print(vciCol, 'VCI Image Collection');

// Aggregate Monthly VCI to Yearly VCI
var yearlyVCI = ee.List.sequence(startYear, endYear).map(function(year) {
  var filtered = vciCol.filter(ee.Filter.calendarRange(year, year, 'year'));
  var yearlyMean = filtered.mean().set('year', year).set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis());
  return yearlyMean;
});

var yearlyVCICol = ee.ImageCollection.fromImages(yearlyVCI);
print(yearlyVCICol, 'Yearly VCI Image Collection');

// ------------------------------------------------- calculation of TCI -----------------------------------------------

var startYear = 2012;
var endYear = 2022;
var startDate = ee.Date.fromYMD(startYear, 1, 1);
var endDate = ee.Date.fromYMD(endYear, 12, 31);
var filtered = modisLST.filter(ee.Filter.date(startDate, endDate));

// Apply QA Mask to select only the highest quality pixels
var bitwiseExtract = function(input, fromBit, toBit) {
  var maskSize = ee.Number(1).add(toBit).subtract(fromBit);
  var mask = ee.Number(1).leftShift(maskSize).subtract(1);
  return input.rightShift(fromBit).bitwiseAnd(mask);
};

var applyQaMask = function(image) {
  var lstDay = image.select('LST_Day_1km');
  var qcDay = image.select('QC_Day');
  var qaMask = bitwiseExtract(qcDay, 0, 1).eq(0); // Highest quality
  var dataQualityMask = bitwiseExtract(qcDay, 2, 3).eq(0);
  var lstErrorMask = bitwiseExtract(qcDay, 6, 7).eq(0);
  var mask = qaMask.and(dataQualityMask).and(lstErrorMask);
  return lstDay.updateMask(mask);
};

var maskedCol = filtered.map(applyQaMask);
var lstCol = maskedCol.select('LST_Day_1km');

// LST scaling
var scaleLST = function(image) {
  var scaled = image.multiply(0.02).subtract(273.15);
  return scaled.copyProperties(image, ['system:index', 'system:time_start']);
};

var lstScaled = lstCol.map(scaleLST);

// Calculate global min and max LST
var globalMin = lstScaled.min();
var globalMax = lstScaled.max();

// Calculate TCI for each image
var calculateTCI = function(image) {
  var tci = image.expression(
    '100 * (max - lst) / (max - min)', {
      'lst': image,
      'max': globalMax,
      'min': globalMin
    }
  ).rename('TCI');
  return tci.copyProperties(image, ['system:time_start']);
};
var tciCol = lstScaled.map(calculateTCI);

var yearlyTCI = ee.List.sequence(startYear, endYear).map(function(year) {
  var filtered = tciCol.filter(ee.Filter.calendarRange(year, year, 'year'));
  var yearlyMean = filtered.mean().set('year', year).set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis());
  return yearlyMean;
});

var yearlyTCICol = ee.ImageCollection.fromImages(yearlyTCI);

// ---------------------------------------------------- calculation of VHI --------------------------------------------

var yearly_merge_vci_tci = yearlyVCICol.map(function(img){
  var sts = img.get('system:time_start')
  var TCI_Band = yearlyTCICol.filter(ee.Filter.eq('system:time_start',sts)).first()
  var VCI_Band = img.select('VCI')
  return VCI_Band.addBands(TCI_Band).copyProperties(img,['system:time_start','year'])
})

print(yearly_merge_vci_tci,'Yearly Merge VCI & TCI')


// Calculate VHI for each image in the collections
var calculateVHI = function(image) {
  var vci = image.select('VCI');
  var tci = image.select('TCI');
  var alpha = 0.5;
  var vhi = vci.expression(
    'alpha * VCI + (1 - alpha) * TCI', {
      'VCI': vci,
      'TCI': tci,
      'alpha': alpha
    }
  ).rename('VHI');
  return vhi.copyProperties(image, ['system:time_start']);
};

var vhiCol = yearly_merge_vci_tci.map(calculateVHI);
print (vhiCol, 'Yearly VhiCol')

// Vegetation Condition Classification  source (Ejaz et al., 2023) 

// | VCI/TCI/VHI          | Condition         |class
// | 0-10                 | EXTREME drought   |1
// | 10-20                | SEVERE drought    |2
// | 20-30                | Moderate drought  |3
// | 30-40                | mild drought      |4
// | more than 40         | no drought        |5


// Function to classify VCI values
var classifyVHI = function(image) {
  var classified = image.select('VHI').expression(
    '(b("VHI") > 40) ? 5 : ' +
    '(b("VHI") > 30) ? 4 : ' +
    '(b("VHI") > 20) ? 3 : ' +
    '(b("VHI") > 10) ? 2 : 1'
  ).rename('drought_class');
  
  return image.addBands(classified);
};

// Apply classification to yearly VCI collection
var classifiedYearlyVHI = vhiCol.map(classifyVHI);
print (classifiedYearlyVHI, 'classifiedYearlyVHI')
// Define a color palette for visualization
var droughtPalette = [
  'FF0000', // Extreme drought (class 1) - Red
  'FFA500', // Severe drought (class 2) - Orange
  'FFFF00', // Moderate drought (class 3) - Yellow
  '98FB98', // Mild drought (class 4) - Pale green
  '006400'  // No drought (class 5) - Dark green
];

// Visualization parameters
var visParams = {
  min: 1,
  max: 5,
  palette: droughtPalette
};

// Function to add layer for a specific year
var addLayerForYear = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year, 12, 31);
  
  var yearImage = classifiedYearlyVHI
    .filter(ee.Filter.date(startDate, endDate))
    .first();
  
  Map.addLayer(yearImage.select('drought_class').clip(aoi), visParams, 'Drought Classification ' + year);
};

// Add layers for each year
for (var year = startYear; year <= endYear; year++) {
  addLayerForYear(year);
}
// Function to export an image for a specific year
var exportImageForYear = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year, 12, 31);
  
  var yearImage = classifiedYearlyVHI
    .filter(ee.Filter.date(startDate, endDate))
    .first();
  
  // Export the image
  Export.image.toDrive({
    image: yearImage.select('drought_class'),
    description: 'Drought_Classification_VHI-' + year,
    folder: 'GEE_Drought_Classification', 
    scale: 1000, // 
    region: aoi, // 
    maxPixels: 1e13
  });
};

// Add layers and start exports for each year
for (var year = startYear; year <= endYear; year++) {
  // addLayerForYear(year);
  exportImageForYear(year);
}

