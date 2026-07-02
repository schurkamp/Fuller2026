# Data for the Study "Improving marsh bird monitoring in Great Lakes wetlands with bioacoustics and BirdNET"

This repository includes the data for Fuller et al. 2026, a study currently in review comparing bioacoustic techniques with in-person surveys to detect secretive marsh birds in the Great Lakes region. A DOI for the study will be linked here upon its publication.

## Metadata: Fuller2026_cleaned_data.csv
- **Route, Point, ARU**: The site, point count location within the site, and name for the recorder of a given point
- **Geography, Lat, Long**: Coordinates of the site along with a categorical variable for its relative location within the Great Lakes, which dictated the survey window
- **Year, Date, StartTime**: Year, date, and time of detection-- either by ARU or within a point count
- **Protocol**:  whether the detection is from ARUs or point counts
- **CommonName, BirdCd**: the species of bird detected and its alpha code
- **Survey**: The survey period corresponding to an observation. Value is either 1, 2, or 3 depending on site and date
- **Visit**: an index for consecutive visits. There are 3 for PCs and 15 for ARUs
- **IndivCount, Distance**: When *Protocol = ARU*, 0 in IndivCount is the species was not found at that time and 1 refers to a positive detection. When *Protocol - PC*, 0 in IndivCount means the species was not found and any value above zero refers to the amount of individuals found during the point count. In that case, Distance is an estimate in meters of how far away the individual was from the point.
- **Added**: Whether the data was actually recorded (value is "No") or added during zero-filling for occupancy analysis (value is "Yes")
- **Confidence**: the score BirdNET-Analyzer gave to its confidence in being correct in detecting the species. Confidence values are only present in rows where *Protocol = ARU*

## Metadata: historical_pc_data.csv
- **TransectName, PointName**: The site and point count location within the site
- **Date, Year, StartTime**: Date, year, and the time the survey playback began
- **Visit**: an index for consecutive visits
- **CommonName, BirdCd, ScientificName**: the species and alpha code for the detected species
- **BirdCount**: number of individuals of the detected species found during the point count
- **IndivId**: a unique observation number
- **DuplicateIndiv**: any non-blank cell in this column represents an observation that needs to be dropped, as it represents a duplicate record
