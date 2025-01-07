// Load image collection and study area 
// var geometry: study area boundary
// var chirps :ImageCollection "CHIRPS Pentad: Climate Hazards Center InfraRed Precipitation With Station Data (Version 2.0 Final)"

// List of 10 years
var lpaYears = ee.List.sequence(2010, 2020)
var months = ee.List.sequence(1, 12)
// Map over the years and create a monthly totals collection
var monthlyImages = lpaYears.map(function(year) {
  return months.map(function(month) {
    var filtered = chirps
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .filter(ee.Filter.calendarRange(month, month, 'month'))
    var monthly = filtered.sum().clip(geometry)
    return monthly.set({'month': month, 'year': year})
  })
}).flatten()

// Create monthly collection
var monthlyCol = ee.ImageCollection.fromImages(monthlyImages)

// mean precipiataion of each month
var longTermMeans = months.map(function(month) {
    var filtered = monthlyCol.filter(ee.Filter.eq('month', month))
    var monthlyMean = filtered.mean().clip(geometry)
    return monthlyMean.set('month', month)
})
var monthlyRainfall = ee.ImageCollection.fromImages(longTermMeans)
// Filter for 2020
var filtered = chirps
  .filter(ee.Filter.date('2020-01-01', '2020-12-31'))
  .filter(ee.Filter.bounds(geometry)) 

// Calculate monthly average rainfall
var monthlyTotals = months
  .map(function(month) {
    return filtered
      .filter(ee.Filter.calendarRange(month, month, 'month'))
        .sum()
        .set('month', month);
});
var observedRainfall = ee.ImageCollection.fromImages(monthlyTotals)
print(observedRainfall)
var palette = ['white', 'blue']
var visParams = {
  min:0,
  max: 2500,
  palette: palette
}
Map.addLayer(monthlyRainfall.sum().clip(geometry),visParams, 'Long Term')
Map.addLayer(observedRainfall.sum().clip(geometry), visParams , 'Current')

// deviation from longterm mean
var deviation = months.map(function(month) {
  var longTermMean = monthlyRainfall
    .filter(ee.Filter.eq('month', month)).first()
  var monthlyObserved = observedRainfall
    .filter(ee.Filter.eq('month', month)).first()
  var deviation = (monthlyObserved.subtract(longTermMean)
    .divide(longTermMean)).multiply(100)
    .set('month', month)
  return deviation
})
// Charting
var chart = ui.Chart.image.series({
  imageCollection: deviation, 
  region: geometry.geometry(), 
  reducer: ee.Reducer.mean(), 
  scale: 5000,
  xProperty: 'month'
}).setOptions({
      interpolateNulls: true,
      lineWidth: 1,
      pointSize: 3,
      title: 'Rainfall Deviation from Long-term Mean',
      vAxis: {title: 'Deviation %'},
      hAxis: {title: 'Month', gridlines: {count: 12}}
});
print(chart);

//  some sections of the scripts in this repository are inpired by open-source learning materials from www.spatialthoughts.com

