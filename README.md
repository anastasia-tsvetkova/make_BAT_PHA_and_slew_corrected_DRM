# make_BAT_PHA_and_slew_corrected_DRM
A class in python and a jupyter notebook providing an example of usage.
The class includes methods to make pha- and rsp-files from Swift/BAT tte files.

It is supposed that you already have a 'standard' folder with BAT data downloaded.

You probably need to run 

> export HEADASPROMPT="/dev/null"

and

> heainit

which is alias to 

> .$HEADAS/headas-init.sh

in your terminal before starting jupyter notebook.

dT0 and ToF are necessary only when obtaining the data relative to the Konus-Wind trigger times.
