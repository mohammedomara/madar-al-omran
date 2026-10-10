// ============================================================================
// MADAR AL-OMRAN | مدار العمران
// V3.2.0 URBAN-GROWTH CORRECTION + ENMAP HYPERSPECTRAL VALIDATION
//
// Why this version exists
// - The published dashboard reported only 0.107 km2 of new urban growth.
// - That estimate was produced by a single hard classification threshold,
//   strict Dynamic World labels and a 3-pixel connected-component filter.
// - V3.1 does NOT choose a larger number by hand. It reports three defensible
//   scenarios, trains the model on two time windows, and adds an independent
//   GHSL plausibility check.
//
// Interpretation rules
// 1) The standard temporally persistent mask is the provisional dashboard mask.
// 2) Strict-to-relaxed persistence is the reported uncertainty range.
// 3) GHSL 2015-2025 is a 100 m modelled/interpolated-or-extrapolated check.
//    It is not a substitute for visual validation or local authoritative data.
// 4) Priority >=1 C is a planning/watch layer; >=2 C is a high-priority layer.
// 5) All output masks are explicitly clipped to the project AOI.
// ============================================================================

// ------------------------------ 0. SETTINGS ---------------------------------

var CONFIG = {
  aoi: ee.Geometry.Rectangle([54.65, 24.31, 54.81, 24.48], null, false),

  hotMonthStart: 5,
  hotMonthEnd: 9,

  historicalStart: '2013-01-01',
  historicalEnd:   '2015-01-01',

  earlyTrainingStart: '2016-01-01',
  earlyTrainingEnd:   '2018-01-01',
  recentTrainingStart: '2023-01-01',
  recentTrainingEnd:   '2025-01-01',

  currentStart: '2025-01-01',
  currentEnd:   '2026-10-01',

  // Automatic reference labels. The relative margin prevents ambiguous
  // built/bare pixels from entering the training data.
  dwBuiltLabelThreshold: 0.45,
  dwBareLabelThreshold: 0.55,
  dwClassMargin: 0.10,

  samplesPerClassPerPeriod: 1200,
  trainFraction: 0.70,
  randomSeed: 813,
  trees: 350,

  // RF probability scenarios. Reporting all three prevents threshold picking.
  conservative: {
    currentMin: 0.70,
    historicalMax: 0.30,
    probabilityChangeMin: 0.30,
    minimumConnectedPixels: 2
  },
  balanced: {
    currentMin: 0.55,
    historicalMax: 0.45,
    probabilityChangeMin: 0.10,
    minimumConnectedPixels: 2
  },
  sensitive: {
    currentMin: 0.40,
    historicalMax: 0.60,
    probabilityChangeMin: 0.05,
    minimumConnectedPixels: 1
  },

  watchPriorityDeltaC: 1.0,
  highPriorityDeltaC: 2.0,

  scale: 30,
  maxPixels: 1e9,
  tileScale: 4,

  // Production mode: keep diagnostic Console requests off so the public app
  // spends its quota on the synchronized dashboard metrics and map tiles.
  debug: false,
  checkEmitCoverage: false,

  // EnMAP L2A surface-reflectance cube acquired on 2025-09-27.
  // The 224 bands are stored as scaled Int16 reflectance (scale 0.0001).
  // EnMAP is used only as independent optical/spectral validation; all LST
  // results remain derived from the Landsat thermal band.
  useHyperspectral: true,
  hyperspectralAsset:
    'projects/tidal-solstice-437700-a3/assets/enmap_shamkha_20250927'
};

var AOI = CONFIG.aoi;

// ------------------------- 1. LANDSAT PREPROCESSING --------------------------

function prepLandsat8(image) {
  var qa = image.select('QA_PIXEL');

  var clear = qa.bitwiseAnd(1 << 0).eq(0)
    .and(qa.bitwiseAnd(1 << 1).eq(0))
    .and(qa.bitwiseAnd(1 << 2).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0));

  var unsaturated = image.select('QA_RADSAT').eq(0);

  var sr = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])
    .multiply(0.0000275)
    .add(-0.2)
    .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);

  var lst = image.select('ST_B10')
    .multiply(0.00341802)
    .add(149.0)
    .subtract(273.15)
    .rename('LST_C');

  var ndvi = sr.normalizedDifference(['nir', 'red']).rename('NDVI');
  var ndbi = sr.normalizedDifference(['swir1', 'nir']).rename('NDBI');
  var mndwi = sr.normalizedDifference(['green', 'swir1']).rename('MNDWI');

  var bsi = sr.expression(
    '((S1 + R) - (N + B)) / ((S1 + R) + (N + B))', {
      S1: sr.select('swir1'),
      R: sr.select('red'),
      N: sr.select('nir'),
      B: sr.select('blue')
    }).rename('BSI');

  return sr.addBands([lst, ndvi, ndbi, mndwi, bsi])
    .updateMask(clear)
    .updateMask(unsaturated)
    .copyProperties(image, ['system:time_start']);
}

function landsatCollection(start, end) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(AOI)
    .filterDate(start, end)
    .filter(ee.Filter.calendarRange(
      CONFIG.hotMonthStart,
      CONFIG.hotMonthEnd,
      'month'
    ))
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .map(prepLandsat8);
}

function seasonalComposite(start, end, name) {
  var collection = landsatCollection(start, end);
  return collection.median().clip(AOI)
    .set('period_name', name)
    .set('scene_count', collection.size());
}

var historical = seasonalComposite(
  CONFIG.historicalStart, CONFIG.historicalEnd, '2013-2014'
);
var earlyTraining = seasonalComposite(
  CONFIG.earlyTrainingStart, CONFIG.earlyTrainingEnd, '2016-2017'
);
var recentTraining = seasonalComposite(
  CONFIG.recentTrainingStart, CONFIG.recentTrainingEnd, '2023-2024'
);
var current = seasonalComposite(
  CONFIG.currentStart, CONFIG.currentEnd, '2025-2026'
);

// Annual composites used only to test that detected change is persistent and
// is not caused by one anomalous seasonal composite.
var historical2013 = seasonalComposite('2013-01-01', '2014-01-01', '2013');
var historical2014 = seasonalComposite('2014-01-01', '2015-01-01', '2014');
var current2025 = seasonalComposite('2025-01-01', '2026-01-01', '2025');
var current2026 = seasonalComposite('2026-01-01', '2026-10-01', '2026');

if (CONFIG.debug) {
  print('=== INPUT SCENE COUNTS ===');
  print('Historical Landsat scenes', historical.get('scene_count'));
  print('Early-training Landsat scenes', earlyTraining.get('scene_count'));
  print('Recent-training Landsat scenes', recentTraining.get('scene_count'));
  print('Current Landsat scenes', current.get('scene_count'));
  print('2013 annual Landsat scenes', historical2013.get('scene_count'));
  print('2014 annual Landsat scenes', historical2014.get('scene_count'));
  print('2025 annual Landsat scenes', current2025.get('scene_count'));
  print('2026 annual Landsat scenes', current2026.get('scene_count'));
}

// --------------------- 2. TWO-WINDOW REFERENCE LABELS -----------------------

function dynamicWorldProbabilities(start, end) {
  return ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
    .filterBounds(AOI)
    .filterDate(start, end)
    .filter(ee.Filter.calendarRange(
      CONFIG.hotMonthStart,
      CONFIG.hotMonthEnd,
      'month'
    ))
    .select(['built', 'bare'])
    .median()
    .clip(AOI);
}

function dynamicWorldLabel(start, end) {
  var probabilities = dynamicWorldProbabilities(start, end);
  var built = probabilities.select('built');
  var bare = probabilities.select('bare');

  var confidentBuilt = built.gte(CONFIG.dwBuiltLabelThreshold)
    .and(built.subtract(bare).gte(CONFIG.dwClassMargin));
  var confidentBare = bare.gte(CONFIG.dwBareLabelThreshold)
    .and(bare.subtract(built).gte(CONFIG.dwClassMargin));

  return ee.Image(0).rename('label')
    .where(confidentBuilt, 1)
    .updateMask(confidentBuilt.or(confidentBare))
    .toByte()
    .clip(AOI);
}

var FEATURES = [
  'blue', 'green', 'red', 'nir', 'swir1', 'swir2',
  'NDVI', 'NDBI', 'MNDWI', 'BSI'
];

function periodSamples(composite, start, end, seedOffset) {
  var stack = composite.select(FEATURES)
    .addBands(dynamicWorldLabel(start, end));

  return stack.stratifiedSample({
    numPoints: CONFIG.samplesPerClassPerPeriod,
    classBand: 'label',
    region: AOI,
    scale: CONFIG.scale,
    seed: CONFIG.randomSeed + seedOffset,
    geometries: true,
    tileScale: CONFIG.tileScale
  });
}

var earlySamples = periodSamples(
  earlyTraining,
  CONFIG.earlyTrainingStart,
  CONFIG.earlyTrainingEnd,
  0
);
var recentSamples = periodSamples(
  recentTraining,
  CONFIG.recentTrainingStart,
  CONFIG.recentTrainingEnd,
  17
);

var samples = earlySamples.merge(recentSamples)
  .randomColumn('split', CONFIG.randomSeed);
var trainSamples = samples.filter(ee.Filter.lt('split', CONFIG.trainFraction));
var testSamples = samples.filter(ee.Filter.gte('split', CONFIG.trainFraction));

var hardClassifier = ee.Classifier.smileRandomForest({
  numberOfTrees: CONFIG.trees,
  seed: CONFIG.randomSeed
}).train({
  features: trainSamples,
  classProperty: 'label',
  inputProperties: FEATURES
});

var testPredictions = testSamples.classify(hardClassifier);
var confusion = testPredictions.errorMatrix('label', 'classification');

// MULTIPROBABILITY returns [P(class 0), P(class 1)] for the two class labels.
var probabilityClassifier = hardClassifier.setOutputMode('MULTIPROBABILITY');

function builtProbability(composite) {
  return composite.select(FEATURES)
    .classify(probabilityClassifier)
    .arrayFlatten([['p_bare', 'p_built']])
    .select('p_built')
    .rename('built_probability')
    .clip(AOI);
}

var historicalBuiltProbability = builtProbability(historical);
var currentBuiltProbability = builtProbability(current);
var historical2013BuiltProbability = builtProbability(historical2013);
var historical2014BuiltProbability = builtProbability(historical2014);
var current2025BuiltProbability = builtProbability(current2025);
var current2026BuiltProbability = builtProbability(current2026);
var builtProbabilityChange = currentBuiltProbability
  .subtract(historicalBuiltProbability)
  .rename('built_probability_change')
  .clip(AOI);

if (CONFIG.debug) {
  print('=== CLASSIFIER VALIDATION ===');
  print('Early-period samples', earlySamples.size());
  print('Recent-period samples', recentSamples.size());
  print('Total samples', samples.size());
  print('Training samples', trainSamples.size());
  print('Testing samples', testSamples.size());
  print('Confusion matrix', confusion);
  print('Overall accuracy', confusion.accuracy());
  print('Kappa', confusion.kappa());
  print('Producer accuracy', confusion.producersAccuracy());
  print('Consumer accuracy', confusion.consumersAccuracy());
}

// -------------------- 3. THREE URBAN-GROWTH SCENARIOS -----------------------

function cleanConnected(mask, minimumPixels) {
  var binary = mask.unmask(0).clip(AOI);
  var connected = binary.selfMask().connectedPixelCount(256, true);
  return binary
    .and(connected.gte(minimumPixels))
    .selfMask()
    .clip(AOI);
}

function urbanGrowthScenario(settings, name) {
  // Water and dense vegetation are excluded, but no bare-soil or BSI rule is
  // imposed here because those rules removed valid paved surfaces previously.
  var plausibleCurrentSurface = current.select('NDVI').lt(0.40)
    .and(current.select('MNDWI').lt(0.10));

  var candidate = currentBuiltProbability.gte(settings.currentMin)
    .and(historicalBuiltProbability.lte(settings.historicalMax))
    .and(builtProbabilityChange.gte(settings.probabilityChangeMin))
    .and(plausibleCurrentSurface);

  return cleanConnected(candidate, settings.minimumConnectedPixels)
    .rename(name)
    .clip(AOI);
}

var newUrbanConservative = urbanGrowthScenario(
  CONFIG.conservative, 'new_urban_conservative'
);
var newUrbanBalanced = urbanGrowthScenario(
  CONFIG.balanced, 'new_urban_balanced'
);
var newUrbanSensitive = urbanGrowthScenario(
  CONFIG.sensitive, 'new_urban_sensitive'
);

// Require the conservative two-year mask to be absent in BOTH historical
// years and present in BOTH current years. This removes pixels caused by one
// anomalous year and produced a substantially more defensible working result.
var newUrbanPersistentStrict = newUrbanConservative
  .and(historical2013BuiltProbability.lte(0.30))
  .and(historical2014BuiltProbability.lte(0.30))
  .and(current2025BuiltProbability.gte(0.70))
  .and(current2026BuiltProbability.gte(0.70))
  .selfMask()
  .rename('new_urban_persistent_strict')
  .clip(AOI);

var newUrbanPersistentStandard = newUrbanConservative
  .and(historical2013BuiltProbability.lte(0.35))
  .and(historical2014BuiltProbability.lte(0.35))
  .and(current2025BuiltProbability.gte(0.65))
  .and(current2026BuiltProbability.gte(0.65))
  .selfMask()
  .rename('new_urban_persistent_standard')
  .clip(AOI);

var newUrbanPersistentRelaxed = newUrbanConservative
  .and(historical2013BuiltProbability.lte(0.40))
  .and(historical2014BuiltProbability.lte(0.40))
  .and(current2025BuiltProbability.gte(0.60))
  .and(current2026BuiltProbability.gte(0.60))
  .selfMask()
  .rename('new_urban_persistent_relaxed')
  .clip(AOI);

// This is the working dashboard layer. Original conservative/balanced/
// sensitive masks remain available as diagnostics, not headline results.
var newUrban = newUrbanPersistentStandard;

// Direct Dynamic World change is shown as corroborating evidence only because
// Dynamic World does not cover the 2013 start year.
var dwEarly = dynamicWorldProbabilities(
  CONFIG.earlyTrainingStart, CONFIG.earlyTrainingEnd
).select('built');
var dwCurrent = dynamicWorldProbabilities(
  CONFIG.currentStart, CONFIG.currentEnd
).select('built');
var dwGrowthEvidence = dwCurrent.gte(0.30)
  .and(dwCurrent.subtract(dwEarly).gte(0.15))
  .selfMask()
  .rename('dw_growth_evidence')
  .clip(AOI);

// --------------------- 4. INDEPENDENT GHSL PLAUSIBILITY CHECK ---------------

// GHSL units are square metres of built surface per 100 m grid cell.
// The 2025 layer may include temporal interpolation/extrapolation; therefore it
// is a broad independent check, not a calibration target and not ground truth.
var ghsl2015 = ee.Image('JRC/GHSL/P2023A/GHS_BUILT_S/2015')
  .select('built_surface')
  .clip(AOI);
var ghsl2025 = ee.Image('JRC/GHSL/P2023A/GHS_BUILT_S/2025')
  .select('built_surface')
  .clip(AOI);
var ghslBuiltSurfaceIncrease = ghsl2025.subtract(ghsl2015)
  .max(0)
  .rename('ghsl_built_surface_increase_m2')
  .clip(AOI);

// Cells with at least 500 m2 modelled increase are displayed for orientation.
var ghslGrowthCells = ghslBuiltSurfaceIncrease.gte(500)
  .selfMask()
  .rename('ghsl_growth_cells')
  .clip(AOI);

var ghslIncreaseKm2 = ee.Number(ghslBuiltSurfaceIncrease.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: AOI,
  scale: 100,
  maxPixels: CONFIG.maxPixels,
  tileScale: CONFIG.tileScale
}).get('ghsl_built_surface_increase_m2')).divide(1e6);

// --------------------------- 5. THERMAL CHANGE -------------------------------

var lstHistorical = historical.select('LST_C').rename('LST_historical');
var lstCurrent = current.select('LST_C').rename('LST_current');
var heatDelta = lstCurrent.subtract(lstHistorical)
  .rename('delta_LST_C')
  .clip(AOI);

var stableDesert = historicalBuiltProbability.lt(0.25)
  .and(currentBuiltProbability.lt(0.25))
  .and(historical.select('NDVI').lt(0.20))
  .and(current.select('NDVI').lt(0.20))
  .and(historical.select('MNDWI').lt(0.00))
  .and(current.select('MNDWI').lt(0.00))
  .selfMask()
  .rename('stable_desert')
  .clip(AOI);

var priorityWatch = newUrban
  .and(heatDelta.gte(CONFIG.watchPriorityDeltaC))
  .selfMask()
  .rename('priority_watch_ge_1C')
  .clip(AOI);

var priorityHigh = newUrban
  .and(heatDelta.gte(CONFIG.highPriorityDeltaC))
  .selfMask()
  .rename('priority_high_ge_2C')
  .clip(AOI);

// ------------------------------ 6. STATISTICS --------------------------------

function maskedStats(image, mask, labelText) {
  var reduced = image.updateMask(mask).reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.median(), '', true)
      .combine(ee.Reducer.stdDev(), '', true)
      .combine(ee.Reducer.count(), '', true)
      .combine(ee.Reducer.percentile([5, 25, 75, 95]), '', true),
    geometry: AOI,
    scale: CONFIG.scale,
    maxPixels: CONFIG.maxPixels,
    tileScale: CONFIG.tileScale
  });
  if (CONFIG.debug) {
    print(labelText, reduced);
  }
  return reduced;
}

function maskAreaKm2(mask) {
  return ee.Number(ee.Image.pixelArea().updateMask(mask).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: AOI,
    scale: CONFIG.scale,
    maxPixels: CONFIG.maxPixels,
    tileScale: CONFIG.tileScale
  }).get('area')).divide(1e6);
}

function validPixelCount(image, mask) {
  return ee.Number(image.updateMask(mask).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: AOI,
    scale: CONFIG.scale,
    maxPixels: CONFIG.maxPixels,
    tileScale: CONFIG.tileScale
  }).values().get(0));
}

var historicalMean = ee.Number(lstHistorical.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: AOI,
  scale: CONFIG.scale,
  maxPixels: CONFIG.maxPixels,
  tileScale: CONFIG.tileScale
}).get('LST_historical'));

var currentMean = ee.Number(lstCurrent.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: AOI,
  scale: CONFIG.scale,
  maxPixels: CONFIG.maxPixels,
  tileScale: CONFIG.tileScale
}).get('LST_current'));

var conservativeAreaKm2 = maskAreaKm2(newUrbanConservative);
var balancedAreaKm2 = maskAreaKm2(newUrbanBalanced);
var sensitiveAreaKm2 = maskAreaKm2(newUrbanSensitive);
var persistentStrictAreaKm2 = maskAreaKm2(newUrbanPersistentStrict);
var persistentStandardAreaKm2 = maskAreaKm2(newUrbanPersistentStandard);
var persistentRelaxedAreaKm2 = maskAreaKm2(newUrbanPersistentRelaxed);
var stableDesertAreaKm2 = maskAreaKm2(stableDesert);
var priorityWatchAreaKm2 = maskAreaKm2(priorityWatch);
var priorityHighAreaKm2 = maskAreaKm2(priorityHigh);
var ghslFootprint100Km2 = maskAreaKm2(
  ghslBuiltSurfaceIncrease.gte(100).selfMask()
);
var ghslFootprint250Km2 = maskAreaKm2(
  ghslBuiltSurfaceIncrease.gte(250).selfMask()
);
var ghslFootprint500Km2 = maskAreaKm2(
  ghslBuiltSurfaceIncrease.gte(500).selfMask()
);
var ghslFootprint1000Km2 = maskAreaKm2(
  ghslBuiltSurfaceIncrease.gte(1000).selfMask()
);

var urbanStats = maskedStats(
  heatDelta, newUrban,
  'Temporally persistent new-urban heat-change statistics'
);
var desertStats = maskedStats(
  heatDelta, stableDesert, 'Stable-desert heat-change statistics'
);

var urbanPixelCount = validPixelCount(heatDelta, newUrban);
var watchPixelCount = validPixelCount(heatDelta, priorityWatch);
var highPixelCount = validPixelCount(heatDelta, priorityHigh);
var urbanDeltaMean = ee.Number(urbanStats.get('delta_LST_C_mean'));
var desertDeltaMean = ee.Number(desertStats.get('delta_LST_C_mean'));
var additionalUrbanHeat = urbanDeltaMean.subtract(desertDeltaMean);

if (CONFIG.debug) {
  print('=== V3.1 CORRECTED RESULTS ===');
  print('Historical AOI mean LST (C)', historicalMean);
  print('Current AOI mean LST (C)', currentMean);
  print('AOI mean LST change (C)', currentMean.subtract(historicalMean));
  print('New urban CONSERVATIVE area (km2)', conservativeAreaKm2);
  print('New urban BALANCED area (km2)', balancedAreaKm2);
  print('New urban SENSITIVE area (km2)', sensitiveAreaKm2);
  print('Persistent new urban STRICT area (km2)', persistentStrictAreaKm2);
  print('Persistent new urban STANDARD area (km2)', persistentStandardAreaKm2);
  print('Persistent new urban RELAXED area (km2)', persistentRelaxedAreaKm2);
  print('GHSL modelled built-surface increase 2015-2025 (km2)', ghslIncreaseKm2);
  print('GHSL growth-cell footprint >=100 m2/cell (km2)', ghslFootprint100Km2);
  print('GHSL growth-cell footprint >=250 m2/cell (km2)', ghslFootprint250Km2);
  print('GHSL growth-cell footprint >=500 m2/cell (km2)', ghslFootprint500Km2);
  print('GHSL growth-cell footprint >=1000 m2/cell (km2)', ghslFootprint1000Km2);
  print('Persistent-standard new-urban valid LST pixels', urbanPixelCount);
  print('New-urban mean heat change (C)', urbanDeltaMean);
  print('Stable-desert area (km2)', stableDesertAreaKm2);
  print('Stable-desert mean heat change (C)', desertDeltaMean);
  print('Additional heat associated with urbanization (C)', additionalUrbanHeat);
  print('Watch-priority >=1 C area (km2)', priorityWatchAreaKm2);
  print('Watch-priority valid pixels', watchPixelCount);
  print('High-priority >=2 C area (km2)', priorityHighAreaKm2);
  print('High-priority valid pixels', highPixelCount);
}

// Create a dedicated ui.Map instance so the dashboard can be mounted safely
// inside ui.SplitPanel after ui.root is cleared.
var dashboardMap = ui.Map();

// ---------------------- 7. HYPERSPECTRAL VALIDATION -------------------------

if (CONFIG.checkEmitCoverage) {
  var emit = ee.ImageCollection('NASA/EMIT/L2A/RFL')
    .filterBounds(AOI)
    .filterDate(CONFIG.currentStart, CONFIG.currentEnd);
  print('=== HYPERSPECTRAL COVERAGE CHECK ===');
  print('NASA EMIT scenes intersecting AOI', emit.size());
  print('EMIT scenes', emit);
}

if (CONFIG.useHyperspectral && CONFIG.hyperspectralAsset !== '') {
  var enmapRaw = ee.Image(CONFIG.hyperspectralAsset).clip(AOI);
  var enmap = enmapRaw.multiply(0.0001);

  // EnMAP band positions nearest to the standard quicklook wavelengths:
  // blue ~463 nm, green ~553 nm, red ~637 nm and NIR ~865 nm.
  // SWIR1/SWIR2 are near 1.6 and 2.2 micrometres. Numeric selectors are
  // zero-based and avoid depending on uploader-generated band names.
  var enmapBlue = enmap.select([9]).rename('enmap_blue');
  var enmapGreen = enmap.select([28]).rename('enmap_green');
  var enmapRed = enmap.select([42]).rename('enmap_red');
  var enmapNir = enmap.select([74]).rename('enmap_nir');
  var enmapSwir1 = enmap.select([131]).rename('enmap_swir1');
  var enmapSwir2 = enmap.select([177]).rename('enmap_swir2');

  var enmapTrueColour = ee.Image.cat([
    enmapRed, enmapGreen, enmapBlue
  ]);
  var enmapFalseColour = ee.Image.cat([
    enmapNir, enmapRed, enmapGreen
  ]);
  var enmapSwirMaterials = ee.Image.cat([
    enmapSwir2, enmapSwir1, enmapRed
  ]);
  var enmapNdvi = enmapNir.subtract(enmapRed)
    .divide(enmapNir.add(enmapRed))
    .rename('EnMAP_NDVI');
  var enmapNdbi = enmapSwir1.subtract(enmapNir)
    .divide(enmapSwir1.add(enmapNir))
    .rename('EnMAP_NDBI');

  if (CONFIG.debug) {
    print('=== ENMAP HYPERSPECTRAL ASSET ===');
    print('EnMAP band count', enmapRaw.bandNames().size());
    print('EnMAP band names', enmapRaw.bandNames());
    print('EnMAP projection', enmapRaw.select(0).projection());
    print('EnMAP metadata', enmapRaw.toDictionary());
  }
}

// ------------------------------ 8. MAP LAYERS --------------------------------

dashboardMap.setOptions('SATELLITE');
dashboardMap.centerObject(AOI, 11);

var rgbVis = {
  bands: ['red', 'green', 'blue'],
  min: 0.02,
  max: 0.35,
  gamma: 1.15
};
var lstVis = {
  min: 25,
  max: 60,
  palette: ['173b8f', '28a7df', 'f7fcb9', 'fdae61', 'd7191c']
};
var deltaVis = {
  min: -5,
  max: 5,
  palette: ['313695', '74add1', 'ffffbf', 'f46d43', 'a50026']
};
var probabilityVis = {
  min: -0.5,
  max: 0.5,
  palette: ['313695', 'ffffbf', 'a50026']
};

dashboardMap.addLayer(historical, rgbVis, 'Satellite 2013-2014', false);
dashboardMap.addLayer(current, rgbVis, 'Satellite 2025-2026', true);
if (CONFIG.useHyperspectral && CONFIG.hyperspectralAsset !== '') {
  dashboardMap.addLayer(
    enmapTrueColour,
    {min: 0.02, max: 0.35, gamma: 1.15},
    'EnMAP 2025 true colour (HSI)', false
  );
  dashboardMap.addLayer(
    enmapFalseColour,
    {min: 0.02, max: 0.45, gamma: 1.10},
    'EnMAP 2025 false colour (NIR-R-G)', false
  );
  dashboardMap.addLayer(
    enmapSwirMaterials,
    {min: 0.02, max: 0.45, gamma: 1.10},
    'EnMAP 2025 SWIR materials (qualitative)', false
  );
  dashboardMap.addLayer(
    enmapNdvi,
    {min: -0.2, max: 0.7, palette: ['7f3b08', 'f6e8c3', '1a9850']},
    'EnMAP 2025 NDVI (HSI)', false
  );
  dashboardMap.addLayer(
    enmapNdbi,
    {min: -0.4, max: 0.4, palette: ['2166ac', 'f7f7f7', 'b2182b']},
    'EnMAP 2025 NDBI (spectral validation)', false
  );
}
dashboardMap.addLayer(lstHistorical, lstVis, 'Historical LST 2013-2014', false);
dashboardMap.addLayer(lstCurrent, lstVis, 'Current LST 2025-2026', false);
dashboardMap.addLayer(heatDelta, deltaVis, 'Long-term heat change', false);
dashboardMap.addLayer(
  builtProbabilityChange, probabilityVis,
  'Built probability change', false
);
dashboardMap.addLayer(
  newUrbanSensitive, {palette: ['ff8c00']},
  'New urban - sensitive', false
);
dashboardMap.addLayer(
  newUrbanBalanced, {palette: ['fff200']},
  'New urban - balanced', false
);
dashboardMap.addLayer(
  newUrbanConservative, {palette: ['00ff66']},
  'New urban - conservative before persistence', false
);
dashboardMap.addLayer(
  newUrbanPersistentRelaxed, {palette: ['ff8c00']},
  'Persistent growth - relaxed', false
);
dashboardMap.addLayer(
  newUrbanPersistentStandard, {palette: ['00ffff']},
  'Persistent growth - standard (working)', true
);
dashboardMap.addLayer(
  newUrbanPersistentStrict, {palette: ['00ff66']},
  'Persistent growth - strict', false
);
dashboardMap.addLayer(
  dwGrowthEvidence, {palette: ['00e5ff']},
  'Dynamic World growth evidence 2016-2026', false
);
dashboardMap.addLayer(
  ghslGrowthCells, {palette: ['cc66ff']},
  'GHSL modelled growth cells 2015-2025', false
);
dashboardMap.addLayer(stableDesert, {palette: ['d9bf77']}, 'Stable desert', false);
dashboardMap.addLayer(
  priorityWatch, {palette: ['ff9800']},
  'Priority watch >=1 C', true
);
dashboardMap.addLayer(
  priorityHigh, {palette: ['ff0000']},
  'High priority >=2 C', true
);
var aoiOutline = ee.Image().byte().paint({
  featureCollection: ee.FeatureCollection([ee.Feature(AOI)]),
  color: 1,
  width: 2
});
dashboardMap.addLayer(aoiOutline, {palette: ['00ffff']}, 'AOI boundary', true);

// ------------------------------ 9. DASHBOARD UI -------------------------------

var panel = ui.Panel({
  style: {
    width: '410px',
    padding: '14px',
    backgroundColor: '#071126'
  }
});

function label(text, size, color, weight) {
  return ui.Label(text, {
    fontSize: size || '13px',
    color: color || '#dce8ff',
    backgroundColor: '#071126',
    fontWeight: weight || 'normal',
    margin: '4px 0'
  });
}

function metricCard(title) {
  var value = ui.Label('Calculating...', {
    fontSize: '21px',
    color: '#ffffff',
    backgroundColor: '#122447',
    fontWeight: 'bold',
    margin: '2px 0'
  });
  panel.add(ui.Panel([
    label(title, '12px', '#9db7e8'),
    value
  ], null, {
    backgroundColor: '#122447',
    padding: '10px',
    margin: '5px 0'
  }));
  return value;
}

panel.add(label('MADAR AL-OMRAN', '24px', '#ffffff', 'bold'));
panel.add(label(
  'مدار العمران | Temporal + Hyperspectral Validation V3.2',
  '16px', '#77d6ff', 'bold'
));
panel.add(label(
  'Al Shamkha – Madinat Al Riyad – Al Falah, Abu Dhabi',
  '12px', '#b7c8e8'
));

var historicalMeanLabel = metricCard('Historical AOI mean LST');
var currentMeanLabel = metricCard('Current AOI mean LST');
var balancedAreaLabel = metricCard('Persistent new urban area - working estimate');
var sensitivityLabel = metricCard('Temporal-persistence uncertainty range');
var ghslLabel = metricCard('GHSL growth-cell footprint (>=250 m2/cell)');
var urbanPixelsLabel = metricCard('Valid thermal pixels in new urban mask');
var urbanDeltaLabel = metricCard('Heat change in new urban areas');
var desertDeltaLabel = metricCard('Heat change in stable desert');
var additionalLabel = metricCard('Urban change relative to stable desert');
var watchAreaLabel = metricCard('Priority watch area >=1 C');
var highAreaLabel = metricCard('High-priority area >=2 C');
var accuracyLabel = metricCard('Pixel hold-out accuracy (not spatial validation)');

panel.add(label('METHOD AND LIMITS', '15px', '#ffffff', 'bold'));
panel.add(label(
  'Landsat 8 L2, May-September, 30 m. Random Forest trained on independent 2016-2017 and 2023-2024 Dynamic World labels.',
  '11px', '#b7c8e8'
));
panel.add(label(
  'A pixel must be absent in both 2013 and 2014, then present in both 2025 and 2026. Strict-standard-relaxed persistence is reported so the result cannot depend on one hidden threshold.',
  '11px', '#b7c8e8'
));
panel.add(label(
  'GHSL is a 100 m modelled plausibility check; its 2025 value may be interpolated/extrapolated and is not ground truth.',
  '10px', '#b7c8e8'
));
panel.add(label(
  'A priority area is a planning threshold, not a claim of statistical significance. Pixel counts are shown for evidence strength.',
  '10px', '#b7c8e8'
));
panel.add(label(
  'The accuracy figure uses a random pixel hold-out from the training windows; it can be optimistic because neighbouring pixels are spatially related.',
  '10px', '#b7c8e8'
));
panel.add(label('HYPERSPECTRAL VALIDATION', '15px', '#ffffff', 'bold'));
panel.add(label(
  'DLR EnMAP HSI L2A surface reflectance, acquired 27 September 2025: 224 contiguous optical bands at 30 m, with 0% cloud cover over the selected scene.',
  '11px', '#b7c8e8'
));
panel.add(label(
  'EnMAP true-colour, false-colour, NDVI, NDBI and SWIR-material layers provide independent current-period spectral evidence. They do not replace or alter the Landsat-derived surface-temperature results.',
  '10px', '#b7c8e8'
));
panel.add(label(
  'Source: DLR EnMAP L2A product DT0000156574, tile 05. Use is governed by the EnMAP data licence; cite DLR (2025), DOI 10.15489/rlyibn8gjc58.',
  '10px', '#8fa8cf'
));

// Resolve every dashboard statistic in one server request. The previous version
// launched a separate evaluate() call for each card; under load, some calls
// could fail or remain pending even though the underlying pixels were valid.
// A single dictionary keeps the cards synchronized and avoids false
// "No valid pixels" messages.
var dashboardMetrics = ee.Dictionary({
  historicalMean: historicalMean,
  currentMean: currentMean,
  persistentStandardArea: persistentStandardAreaKm2,
  persistentStrictArea: persistentStrictAreaKm2,
  persistentRelaxedArea: persistentRelaxedAreaKm2,
  ghslFootprint250: ghslFootprint250Km2,
  urbanPixelCount: urbanPixelCount,
  urbanDeltaMean: urbanDeltaMean,
  desertDeltaMean: desertDeltaMean,
  additionalUrbanHeat: additionalUrbanHeat,
  priorityWatchArea: priorityWatchAreaKm2,
  priorityHighArea: priorityHighAreaKm2,
  accuracyPercent: confusion.accuracy().multiply(100)
});

var metricWidgets = [
  historicalMeanLabel,
  currentMeanLabel,
  balancedAreaLabel,
  sensitivityLabel,
  ghslLabel,
  urbanPixelsLabel,
  urbanDeltaLabel,
  desertDeltaLabel,
  additionalLabel,
  watchAreaLabel,
  highAreaLabel,
  accuracyLabel
];

function setLoadingState(message) {
  metricWidgets.forEach(function(widget) {
    widget.setValue(message);
  });
}

function formatMetric(value, decimals, suffix) {
  if (value === null || value === undefined || !isFinite(Number(value))) {
    return 'Data unavailable';
  }
  return Number(value).toFixed(decimals) + suffix;
}

dashboardMetrics.evaluate(function(metrics, error) {
  if (error || !metrics) {
    setLoadingState('Data unavailable — reload app');
    return;
  }

  historicalMeanLabel.setValue(
    formatMetric(metrics.historicalMean, 2, ' °C')
  );
  currentMeanLabel.setValue(
    formatMetric(metrics.currentMean, 2, ' °C')
  );
  balancedAreaLabel.setValue(
    formatMetric(metrics.persistentStandardArea, 3, ' km²')
  );

  if (
    metrics.persistentStrictArea === null ||
    metrics.persistentStrictArea === undefined ||
    metrics.persistentRelaxedArea === null ||
    metrics.persistentRelaxedArea === undefined
  ) {
    sensitivityLabel.setValue('Data unavailable');
  } else {
    sensitivityLabel.setValue(
      Number(metrics.persistentStrictArea).toFixed(3) + ' – ' +
      Number(metrics.persistentRelaxedArea).toFixed(3) + ' km²'
    );
  }

  ghslLabel.setValue(
    formatMetric(metrics.ghslFootprint250, 3, ' km²')
  );
  urbanPixelsLabel.setValue(
    formatMetric(metrics.urbanPixelCount, 0, ' pixels')
  );
  urbanDeltaLabel.setValue(
    formatMetric(metrics.urbanDeltaMean, 2, ' °C')
  );
  desertDeltaLabel.setValue(
    formatMetric(metrics.desertDeltaMean, 2, ' °C')
  );
  additionalLabel.setValue(
    formatMetric(metrics.additionalUrbanHeat, 2, ' °C')
  );
  watchAreaLabel.setValue(
    formatMetric(metrics.priorityWatchArea, 3, ' km²')
  );
  highAreaLabel.setValue(
    formatMetric(metrics.priorityHighArea, 3, ' km²')
  );
  accuracyLabel.setValue(
    formatMetric(metrics.accuracyPercent, 1, ' %')
  );
});

ui.root.clear();
ui.root.add(ui.SplitPanel({
  firstPanel: panel,
  secondPanel: dashboardMap,
  orientation: 'horizontal',
  wipe: false,
  style: {stretch: 'both'}
}));

// --------------------------- 10. EXPORTS (OPTIONAL) --------------------------

var validationStack = ee.Image.cat([
  historicalBuiltProbability,
  currentBuiltProbability,
  builtProbabilityChange,
  newUrbanConservative.unmask(0),
  newUrbanBalanced.unmask(0),
  newUrbanSensitive.unmask(0),
  newUrbanPersistentStrict.unmask(0),
  newUrbanPersistentStandard.unmask(0),
  newUrbanPersistentRelaxed.unmask(0),
  heatDelta,
  priorityWatch.unmask(0),
  priorityHigh.unmask(0)
]).clip(AOI);

// Uncomment only after inspecting all three masks against the two satellite
// composites and recording the selected scenario in the PoC methodology.
/*
Export.image.toDrive({
  image: validationStack,
  description: 'Madar_Al_Omran_V3_1_Urban_Growth_Validation',
  folder: 'Madar_Al_Omran',
  region: AOI,
  scale: CONFIG.scale,
  maxPixels: CONFIG.maxPixels
});
*/
