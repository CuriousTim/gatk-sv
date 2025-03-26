###########################
# IMPORT LIBRARIES
###########################

import argparse
import pandas as pd
import numpy as np
import math

###########################
# PARSE INPUT ARGUMENTS
###########################

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('--bed', dest='bed', help='Input BED file')
parser.add_argument('--denovo_outlier_factor', dest='denovo_outlier_factor', help='Factor by which to multiply Q3 + 1.5*IQR to define the threshold to call a sample an outlier')
args = parser.parse_args()

bed_file = args.bed
denovo_outlier_factor = args.denovo_outlier_factor

###########################
# Load data
###########################

bed_original = pd.read_csv(bed_file, sep='\t', index_col=False).replace(np.nan, '', regex=True)

###########################
# ANALYSES
###########################

# Remove outliers
bed = bed_original
sample_counts = bed['sample'].value_counts()
q3, q1 = np.percentile(sample_counts.tolist(), [75, 25])
iqr = q3 - q1
count_threshold = int(math.ceil((1.5 * iqr + q3) * denovo_outlier_factor))
sample_counts2 = sample_counts.rename_axis('sample').reset_index(name='count')
samples_keep = sample_counts2[(sample_counts2['count'] <= count_threshold)]['sample'].tolist()

# Keep samples and outliers in sepparate files
output = bed[(bed['sample'].isin(samples_keep))]
output_outliers = bed[(~bed['sample'].isin(samples_keep))]

bed_original.loc[~(bed_original['sample'].isin(samples_keep)) & bed_original['is_de_novo'], 'filter_flag'] = 'sample_outlier'

###########################
# WRITE OUTPUT
###########################

output.to_csv(path_or_buf='final.denovo.merged.bed', mode='a', index=False, sep='\t', header=True)
output_outliers.to_csv(path_or_buf='final.denovo.merged.outliers.bed', mode='a', index=False, sep='\t', header=True)

bed_original.to_csv(path_or_buf='final.denovo.merged.allSamples.bed', mode='a', index=False, sep='\t', header=True)
