from grid_data_format import *
import sys
# Assumes that last input is the basename for the athena dataset.
# i.e. kh_3d_mhd_hlld_128_beta5000_sub_tanhd.0030
basename = sys.argv[-1]
converter = AthenaConverter(basename)
converter.convert()
