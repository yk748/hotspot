# hotspot

Hotspot identification is achieved through an algorithm that detects periods where statistical features of passive sensing variables and 
stress indicators change, highlighting intervals of interest. We first extend the moving sum (MOSUM) change-point method to capture simultaneous changes 
both within and across series by incorporating two types of changes using Mahalanobis distances. Hotspots are then defined in two ways: either as collections 
of time points where the distances exceed specified thresholds, which are subsequently combined, or as the timings corresponding to pointwise confidence intervals 
for individual change points obtained via bootstrapping. A fully functioning package is under development and is expected to be released soon. 
The corresponding manuscript is Kim, Basu, and Banerjee ([2025](#ref-hotspot))

## References

<div id="ref-hotspot" class="references">

Younghoon Kim, Sumanta Basu, and Samprit Banerjee. 2025. “Identification of Emotionally Stressful Periods Through Tracking Changes in Statistical Features of mHealth Data.” *Preprint*
