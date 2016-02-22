How to use the droplet counter plugins
======================================

write bug reports, feature requests and questions to:
samimoll@googlemail.com



What you need:
--------------
 * a recent(2008) version of ImageJ. If the plugin gives you strange
   errors, your version of ImageJ is probably too old.
 * the "Droplet_Finder.jar" file

Quick Start:
------------
 * Copy all the needed files into your ImageJ/plugins folder
 * Start ImageJ, open your image stack, select "Droplet Finder/Filterstacker"
   - Enter the minimum and maximum feature sizes (size bounds in pixels for your droplets)
   - Enter your Z/X aspect ratio
 * wait until a Window labeled "BP-..." appears (could take some time)
   - In the "Watershed 3D" Box enter "2.0" for all radii, check "Invert", then "OK"
 * In the "Segment Analyzer" Box choose "BP-..." as image, "WS-..." as mask
   - If you run Windows: before you click ok, make sure that the "Filterstacker-..."-image is
     not occluded by any window, otherwise you won't see the preview
 * navigate the stack with the "Slice" slider while choosing optimal area and maximum
   threshold parameters. Click "OK"
 * The Segment Analyzer will output a measurement table



What the Filterstacker actually does
------------------------------------
The "Filterstacker" plugin actually only does a 3D-bandpass on your stack. It will first
3D-blur the input image with a filter size of "maximum feature size" and substract the
blurred image from the original. It will then eliminate small features (like noise) by
blurring again with a small filter the size of "minimum feature size". The "Z/X aspect
ratio" compensates for different lateral and vertical resolutions.

It will then call the "Watershed 3D" plugin that finds all local maxima (white spots)
and grows regions around them, so that each region only contains one maximum. It will
output a new stack where each region is labeled with a unique color. The radii
control how far apart the maxima must be until they are regarded as only one maximum.
Setting this to higher values makes the watershed transform more resistant to noise
while increasing the running time (higher radii than 3 take VERY long) and the
probability that two spots that are close together are combined into one region. Big
radii also make the edges of the regions more fuzzy. Setting all radii to 2 is a good
compromise. If the "Invert" option is unchecked, the plugin will find minima (black
spots) instead of maxima.

Finally the Filterstacker calls the "Segment Analyzer" plugin that takes as input an
image and a mask (=the watershed transform). The Segment Analyzer decides which regions
contain a droplet by thresholding (explained below). It then does a FWHM threshold on
all regions that passed the previous thresholding test and measures the volume, position
and surface area of each droplet.
The region thresholding is very simple at the moment. It only considers the maximal and
summed (over the whole region) brightness (the "maximum threshold" and "area threshold"
parameters). If both values are above the respective thresholds, the region is considered
to contain a droplet.
The Segment Analyzer has a preview. Because of ImageJ design restrictions, the stack
has to be navigated with the "Slice" slider. Red marks particles, the borders between
particles are green. If there is a green line through one of your particles, lower the
"connect threshold".
