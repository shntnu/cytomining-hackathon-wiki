# Cell segmentation and field-of-view correction 
•	Identify large, bright, non-segmentable  objects in nuclear channel
•	Classify objects as “Nuclear_Debris” for elimination
•	Define well mask based on “Nuclear_Debris” for all channels for debris
•	Apply “Nuclear_Debris” well mask and identify primary objects in nuclear channel as “Nuclei”
•	Basic size & intensity thresholding for valid objects
•	Basic shape-based de-clumping
•	Filter nuclear objects based on form factor/eccentricity  (i.e. roundness) and classify as “Valid_Nuclei”
•	Identify secondary objects via propagation from nuclei
•	Identify secondary objects via propagation from nuclei
•	Measure object size & shape for Nuclei/Actin/Tubulin objects
•	Measure object intensity in Actin and Tubulin images for each of the 3 object classes
•	Measure image intensity for all 3 channels
•	Measure image area occupied for each of the three object classes
•	Measure texture in actin source image for actin objects and in tubulin source image for tubulin objects
•	Measure granularity in Actin and Tubulin source images for their respective objects
•	Export data

@csmolnar
@JCHTB
