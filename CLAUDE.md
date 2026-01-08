# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**ijp-jacop-b** is an ImageJ/Fiji plugin for colocalization analysis in microscopy images. It's a revamped version of the original JACoP plugin by F. Cordelières and S. Bolte, developed by the BIOP (BioImaging And Optics Platform) at EPFL. The plugin analyzes multi-channel microscopy images to measure colocalization between fluorescent signals using various statistical methods.

The plugin provides:
- Multiple colocalization metrics (Pearson's, Mander's, Costes, Li, Spearman Rank)
- ROI management for per-cell or per-region analysis
- Z-stack handling (2D vs 3D analysis options)
- Automatic and manual thresholding methods
- Random Costes analysis for statistical validation
- Output visualizations (fluorograms, thresholded masks, randomized images)

## Build System

This project uses Maven with the SciJava parent POM.

### Build Commands

```bash
# Standard build
mvn clean package

# Build with tests
mvn clean verify

# Install to local Maven repository
mvn clean install

# Deploy (requires credentials)
mvn clean deploy
```

### CI/CD

- GitHub Actions runs builds on push/PR to master branch
- CI uses Java 8 (Zulu distribution)
- Build scripts are downloaded from scijava/scijava-scripts repository
- The workflow is defined in `.github/workflows/build.yml`

### Installation for Development

To test the plugin in Fiji during development, uncomment and set the path in `pom.xml`:
```xml
<imagej.app.directory>C:/path/to/Fiji/</imagej.app.directory>
```

Then run `mvn clean install` to copy the built JAR to your Fiji installation.

## Code Architecture

### Package Structure

```
ch.epfl.biop.coloc
├── JACoP_B.java              # Main plugin entry point (implements PlugIn)
└── utils/
    ├── ImageColocalizer.java  # Core colocalization analysis engine
    ├── RandomCostes.java      # Random Costes statistical validation
    ├── Utils.java             # Image processing utilities
    ├── ImgInfo.java          # Image metadata wrapper
    └── Object3D.java         # 3D object handling (legacy)
```

### Key Components

**JACoP_B** (`ch.epfl.biop.coloc.JACoP_B`)
- Main plugin class that implements ImageJ's `PlugIn` interface
- Handles user dialog, parameter collection, and workflow orchestration
- Supports both single image and batch folder processing
- Manages ROI iteration and Z-slice processing options
- Creates output montages combining masks, fluorograms, and results

**ImageColocalizer** (`ch.epfl.biop.coloc.utils.ImageColocalizer`)
- Core analysis engine performing all colocalization calculations
- Key methods:
  - `setup()` - Initializes image pairs and calibration
  - `setThresholds()` - Sets manual or automatic thresholds
  - `Pearson()` - Pearson's correlation coefficient
  - `SpearmanRank()` - Spearman's rank correlation
  - `Overlap()` - Overlap coefficients
  - `MM()` - Mander's coefficients (M1, M2, tM1, tM2)
  - `CostesAutoThr()` - Costes automatic threshold determination
  - `RandomCostes2D()` - Random shuffling statistical validation
  - `CytoFluo()` - Generates fluorogram (scatter plot)
  - `ICA()` - Intensity Correlation Analysis (Li's method)
- Uses Fiji's Coloc2 algorithms via `sc.fiji.coloc` package
- Results are stored in an internal `ResultsTable`

**RandomCostes** (`ch.epfl.biop.coloc.utils.RandomCostes`)
- Implements the Costes randomization approach for statistical significance testing
- Performs block-wise shuffling of one channel while keeping the other fixed
- Compares observed Pearson correlation to distribution from randomized images
- Calculates p-value for colocalization significance

**Utils** (`ch.epfl.biop.coloc.utils.Utils`)
- Image processing utilities for colocalization workflows
- `clearOutside()` - Custom ROI masking (fixes ImageJ's fillOutside shift bug)
- `binarize()` - Creates binary masks from thresholded images
- Display range management utilities

### Data Flow

1. **User Input** → JACoP_B collects parameters via `GenericDialogPlus`
2. **Image Preparation** → Channels extracted, ROIs applied, Z-slices handled
3. **Analysis Setup** → ImageColocalizer initialized with image pairs
4. **Thresholding** → Manual values or automatic methods (using ImageJ's AutoThreshold)
5. **Colocalization Metrics** → Selected analyses performed by ImageColocalizer
6. **Statistical Validation** → Optional RandomCostes shuffling analysis
7. **Output Generation** → Results table + montage images (masks, fluorogram, etc.)

### Integration with Fiji Ecosystem

- Uses `sc.fiji.coloc` (Coloc2) for core algorithms (Pearson's, Spearman's, Costes auto-threshold)
- Leverages ImgLib2 for efficient pixel access (`RandomAccessibleInterval`)
- Integrates with ImageJ's ROI Manager for multi-region analysis
- Uses ImageJ's auto-threshold methods via `ij.plugin.Thresholder`
- Depends on `Multi_Stack_Montage` (BIOP plugin) for output visualization

### Important Implementation Notes

- **Costes 3D limitation**: Costes auto-threshold currently only works in 2D. For 3D stacks, use "Consider Z slices Separately" option
- **ROI shift bug workaround**: Custom `clearOutside()` implementation in Utils.java fixes a (1,1) pixel shift in ImageJ's `fillOutside()` when used programmatically
- **Stack histogram thresholding**: Can threshold each slice independently or use histogram from entire stack
- **Preference storage**: User settings stored via `ij.Prefs` with prefix `"jacop.b."`

## Plugin Registration

The plugin is registered in `src/main/resources/plugins.config`:
```
Plugins>BIOP>Image Analysis, "BIOP JACoP", ch.epfl.biop.coloc.JACoP_B
```

After installation via the PTBIOP update site, users access it via: `Plugins > BIOP > Image Analysis > BIOP JACoP`

## Scientific Context

This plugin implements methods from these key papers:
- Pearson's correlation coefficient (standard statistical measure)
- Mander's coefficients (Manders et al., 1993)
- Costes automatic threshold (Costes et al., 2004)
- Li's Intensity Correlation Analysis (Li et al., 2004)

Users should always include control samples (mono-stained) to validate threshold settings and account for crosstalk/bleedthrough.