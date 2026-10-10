# Madar Al-Omran | مدار العمران

**Satellite-based urban growth and land-surface-temperature monitoring for Abu Dhabi**

Madar Al-Omran is a proof of concept developed for the National Space Academy 813 Hackathon challenge on urban expansion, land-use change and health risk. It compares urban growth and land surface temperature across Al Shamkha, Madinat Al Riyad and Al Falah, then highlights locations that may deserve closer planning and heat-mitigation attention.

## Live dashboard

[Open the Google Earth Engine dashboard](https://tidal-solstice-437700-a3.projects.earthengine.app/view/madar-al-omran-dashboard)

## Problem

Urban expansion replaces open land with buildings, roads and other surfaces that change how solar energy is absorbed, stored and released. The planning need is not simply to measure how much a city has expanded, but to identify where expansion and surface-temperature change overlap, while separating that signal from the wider regional warming seen in stable desert areas.

## Proposed solution

The workflow:

1. Builds cloud-masked Landsat 8 Level-2 seasonal composites for May–September.
2. Compares the 2013–2014 historical period with the 2025–2026 current period.
3. Trains a Random Forest urban classifier using independent 2016–2017 and 2023–2024 Dynamic World reference windows.
4. Requires detected change to be temporally persistent: absent in both 2013 and 2014, then present in both 2025 and 2026.
5. Reports conservative, standard and sensitive scenarios instead of relying on one hidden threshold.
6. compares land surface temperature change in new urban areas with stable desert.
7. Uses GHSL as an independent plausibility check.
8. Adds EnMAP hyperspectral layers as independent current-period spectral validation.

## Area of interest

- Al Shamkha
- Madinat Al Riyad
- Al Falah
- Abu Dhabi, United Arab Emirates
- Approximate bounding box: 54.65–54.81°E, 24.31–24.48°N

## Latest validated outputs

| Metric | Result |
|---|---:|
| Historical AOI mean LST | 51.38 °C |
| Current AOI mean LST | 53.11 °C |
| Persistent new urban area — working estimate | 26.202 km² |
| Temporal-persistence uncertainty range | 21.537–28.925 km² |
| GHSL growth-cell footprint, ≥250 m²/cell | 30.497 km² |
| Valid thermal pixels in new-urban mask | 32,104 |
| Heat change in new urban areas | +1.16 °C |
| Heat change in stable desert | +1.81 °C |
| Urban change relative to stable desert | −0.65 °C |
| Priority watch area, ≥1 °C | 14.287 km² |
| High-priority area, ≥2 °C | 7.145 km² |
| Random pixel hold-out accuracy | 98.9% |

### Scientific interpretation

The current run detects substantial persistent urban expansion and localized areas of increasing surface temperature. However, stable desert warmed more than the newly urbanized mask during the comparison period. Therefore, this PoC does **not** claim that urbanization caused additional AOI-wide warming beyond the desert background. The priority layers are planning-screening thresholds, not claims of statistical significance or health outcomes.

## Hyperspectral validation

The final code integrates:

- DLR EnMAP Level-2A surface reflectance
- Acquisition: 27 September 2025 at 07:37:39 UTC
- Product: DT0000156574, tile 05
- 224 contiguous optical bands
- 30 m spatial resolution
- 0% reported scene cloud cover
- DOI: [10.15489/rlyibn8gjc58](https://doi.org/10.15489/rlyibn8gjc58)

EnMAP true-colour, false-colour, NDVI, NDBI and qualitative SWIR-material layers provide independent spectral evidence for the current period. They do not replace or alter the Landsat-derived thermal results.

**Transparency note:** this build uses DLR EnMAP data. It does not claim that the hyperspectral scene was acquired by satellite 813.

## Data sources and use conditions

| Dataset | Purpose | Provider / conditions |
|---|---|---|
| Landsat 8 Collection 2 Tier 1 Level 2 | Surface reflectance, indices and thermal LST | USGS/NASA; subject to USGS Landsat data policy |
| Dynamic World V1 | Reference labels for classifier training and corroborating change evidence | Google and World Resources Institute; subject to provider terms |
| GHSL P2023A built surface | Independent 2015–2025 plausibility check | European Commission Joint Research Centre; cite the GHSL product |
| EnMAP L2A | Hyperspectral validation | DLR; governed by the EnMAP data licence and citation requirements |

The repository does not redistribute the EnMAP raster. Users reproducing the workflow must use an authorised Earth Engine asset and comply with the original data-provider licences.

## Repository contents

- `madar_al_omran_v3_2.js` — complete Google Earth Engine dashboard and analysis code.
- `README.md` — methodology, outputs, sources, licences and limitations.

## Running the code

1. Open the [Google Earth Engine Code Editor](https://code.earthengine.google.com/).
2. Copy the contents of `madar_al_omran_v3_2.js` into a new script.
3. Confirm access to the configured EnMAP Earth Engine asset, or replace its asset path with an authorised copy.
4. Run the script and review the map layers, metrics and console diagnostics.
5. For deployment, update the existing Earth Engine App rather than creating a second dashboard URL.

## Limitations

- Landsat and EnMAP are 30 m products; small features and mixed pixels may be unresolved.
- The 2026 composite is partial through 1 October 2026.
- GHSL is a 100 m modelled plausibility check and is not local ground truth.
- Dynamic World does not cover the 2013 start year and is used through independent later training windows.
- Random pixel hold-out accuracy may be optimistic because neighbouring pixels are spatially related.
- Final planning decisions require field validation, local authoritative land-use data and additional statistical analysis.

## Team

**Madar Al-Omran | مدار العمران**

Hackathon team: Mohammed, Youssef and Mariam.
