from GCMtools.GCMT import GCMT

# Construct a GCMtools object
gcmt = GCMT()

# Read in a GCM dataset...
# 1. directly from MITgcm output, or ...
gcmt.read_raw('MITgcm', mydatapath) # TODO
# 2. from previously reduced data
gcmt.read_reduced('GCMtools/example_dataset.nc', tag='HD209458b')

# Have a first quick look at the data
gcmt.plot_overview(tag='HD209458b') # TODO

# Access all models
all_models = gcmt.models
all_modes = gcmt.get_models()
# or a specific one using the tag
ds = gcmt.models['HD209458b']
ds = gcmt.get_models('HD209458b')
