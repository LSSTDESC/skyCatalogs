import yaml
from skycatalogs.utils.config_utils import get_file_metadata
import sys

# args are path to file to be read and text (yaml) file to dump it to

if len(sys.argv) == 1:
    print('to invoke: \n   python dump_meta.py path-to-input-parquet [path-to-output-text]')
    print('If no output filepath is supplied, output will go to std out')
    exit(0)
data = get_file_metadata(sys.argv[1])

if len(sys.argv) > 2:
    with open(sys.argv[2], mode='w') as f:
        yaml.dump(data, f)
else:
    print(yaml.dump(data))
